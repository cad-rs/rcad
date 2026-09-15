// OCCT AppCont (TKGeomBase/AppCont) — 1:1 Rust translation:
// - AppCont_Function.hxx L27-80            (the continuous-function interface)
// - AppCont_LeastSquare.hxx L27-68         (the class declaration)
// - AppCont_LeastSquare.cxx L30-614        (the bodies)
// - AppCont_ContMatrices.pxx L17-181       (the Bernstein matrix providers)
//
// The package is consumed by Approx_FitAndDivide (the Approx_ComputeCLine.gxx
// instantiation with MultiLine = AppCont_Function,
// Approx_FitAndDivide_0.cxx L21-25) — see approx_fit_and_divide.rs.
//
// Architecture difference (the ContMatrices tables only): OCCT stores the
// MMatrix / InvMMatrix / IBPMatrix / IBTMatrix / VBernstein payloads as
// precomputed `static const double` tables (~124k generated data lines across
// AppCont_ContMatrices_{BB,InvM,IBP,IBT,VB}.pxx — too large for the rcad
// 2000-line file rule).  rcad computes the same matrices at runtime from
// their exact closed forms:
//   M(i,j)      = C(n,i)C(n,j) / ((2n+1) C(2n,i+j))   (Bernstein mass matrix,
//                 n = classe-1; exact integer numerators < 2^53)
//   InvM        = M^-1                                (math_Matrix::Inverse)
//   IBP / IBT   = (M_p^T M_p)^-1 / (M_t^T M_t)^-1     (the same runtime
//                 fallback AppCont_LeastSquare.cxx L486-499 uses for
//                 out-of-table classes)
//   VB(i,j)     = B_i(u_j), u_j = (1 - GaussPoint(j))/2  (row-major, the
//                 ordering verified against the OCCT table values)
// The unit tests pin the computed values to literal values extracted from
// the OCCT .pxx tables (classe 9 / nbpoints 17 — the Approx_FitAndDivide
// degree-8 path).

use glam::{DVec2, DVec3};
use rcad_kernel::core::precision::{CONFUSION, REAL_LAST};
use rcad_kernel::math::gauss_points::{gauss_points, gauss_weights};
use rcad_kernel::math::math_matrix::{Matrix, Vector as RVector};
use rcad_kernel::math::VecD;

use super::approx_int::{AppParConstraint, MultiCurve, MultiPoint};

/// OCCT Standard::Epsilon(theValue) (Standard_Real.hxx L242-250) — the
/// distance to the nearest neighbouring double in the direction of infinity
/// with the same sign as theValue.
fn standard_epsilon(the_value: f64) -> f64 {
    if the_value >= 0.0 {
        {
            let next = if the_value == REAL_LAST {
                REAL_LAST
            } else {
                f64::from_bits(the_value.to_bits() + 1)
            };
            next - the_value
        }
    } else {
        the_value - f64::from_bits(the_value.to_bits() - 1)
    }
}

/// OCCT Standard::Epsilon(1.) — the constant used by
/// AppCont_LeastSquare::FixSingleBorderPoint (the .cxx L44).
const EPSILON_1: f64 = f64::EPSILON;

// ---------------------------------------------------------------------------
// OCCT AppCont_Function (AppCont_Function.hxx L27-80).
// ---------------------------------------------------------------------------

/// OCCT AppCont_Function — class describing a continuous 3d and/or 2d
/// function f(u); provided by the user to the approximation algorithm.
///
/// The OCCT virtual dispatch maps to a Rust trait; the `myNbPnt`/`myNbPnt2d`
/// protected fields are carried by the required `my_nb_pnt`/`my_nb_pnt2d`
/// reads.
pub trait AppContFunction {
    /// OCCT: the protected myNbPnt field.
    fn my_nb_pnt(&self) -> i32;
    /// OCCT: the protected myNbPnt2d field.
    fn my_nb_pnt2d(&self) -> i32;

    /// OCCT GetNumberOfPoints(theNbPnt, theNbPnt2d) (hxx L37-41).
    fn get_number_of_points(&self) -> (i32, i32) {
        (self.my_nb_pnt(), self.my_nb_pnt2d())
    }

    /// OCCT GetNbOf3dPoints() (hxx L44).
    fn get_nb_of_3d_points(&self) -> i32 {
        self.my_nb_pnt()
    }

    /// OCCT GetNbOf2dPoints() (hxx L47).
    fn get_nb_of_2d_points(&self) -> i32 {
        self.my_nb_pnt2d()
    }

    /// OCCT FirstParameter() (hxx L53).
    fn first_parameter(&self) -> f64;

    /// OCCT LastParameter() (hxx L56).
    fn last_parameter(&self) -> f64;

    /// OCCT Value(theU, thePnt2d, thePnt) (hxx L59-61).
    fn value(&self, the_u: f64, the_pnt2d: &mut [DVec2], the_pnt: &mut [DVec3]) -> bool;

    /// OCCT D1(theU, theVec2d, theVec) (hxx L64-66).
    fn d1(&self, the_u: f64, the_vec2d: &mut [DVec2], the_vec: &mut [DVec3]) -> bool;

    /// OCCT PeriodInformation(theDimIdx, IsPeriodic, thePeriod) (hxx L71-75)
    /// — the base implementation reports no periodicity.
    fn period_information(&self, _the_dim_idx: i32) -> (bool, f64) {
        (false, 0.0)
    }
}

// ---------------------------------------------------------------------------
// OCCT AppCont_ContMatrices (AppCont_ContMatrices.pxx L17-181).
// ---------------------------------------------------------------------------

/// OCCT THE_MAX_CLASSE (the .pxx L30).
const THE_MAX_CLASSE: i32 = 26;
/// OCCT THE_MAX_MMATRIX_CLASSE (the .pxx L31).
const THE_MAX_MMATRIX_CLASSE: i32 = 24;
/// OCCT THE_MAX_NBPOINTS (the .pxx L32).
const THE_MAX_NBPOINTS: i32 = 24;

/// OCCT AppCont_ContMatrices::MMatrixOffset (the .pxx L39-48) — the packed
/// offset form of the (generated) table layout; rcad keeps it as the index
/// sanity anchor of the runtime builders.
#[allow(dead_code)]
fn m_matrix_offset(the_classe: i32) -> i32 {
    let mut an_offset = 0;
    for i in 2..the_classe {
        an_offset += i * i;
    }
    an_offset
}

/// C(n, k) computed exactly in f64: every intermediate value is an integer
/// below 2^53 for the n <= 46 used by the table range (classe <= 26).
fn binomial(n: i32, k: i32) -> f64 {
    if k < 0 || k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut result = 1.0f64;
    for i in 1..=k {
        // result * (n - k + i) / i stays integral at every step.
        result = result * (n - k + i) as f64 / i as f64;
    }
    result
}

/// The Bernstein mass matrix M(i,j) = int_0^1 B_i(t) B_j(t) dt for the
/// degree n = classe-1 Bernstein basis (the exact rational closed form of
/// the OCCT THE_BBMATRIX_DATA table).
fn bernstein_mass_matrix(classe: i32) -> Matrix {
    let n = classe - 1;
    let mut m = Matrix::new(1, classe, 1, classe);
    for i in 0..classe {
        for j in 0..classe {
            // C(n,i) C(n,j) / ((2n+1) C(2n, i+j)).
            let num = binomial(n, i) * binomial(n, j);
            let den = (2 * n + 1) as f64 * binomial(2 * n, i + j);
            m.set(i + 1, j + 1, num / den);
        }
    }
    m
}

/// The Bernstein basis B_k(u) of degree n (used by VBernstein).
fn bernstein_value(n: i32, k: i32, u: f64) -> f64 {
    binomial(n, k) * u.powi(k) * (1.0 - u).powi(n - k)
}

/// Copies the (classe x classe) 1-based source into the target matrix
/// through the target's own lower bounds — the OCCT .pxx write form
/// `M(i + M.LowerRow() - 1, j + M.LowerCol() - 1) = aData[k++]`.
fn fill_with_bounds(target: &mut Matrix, source: &Matrix) {
    let lr = target.lower_row;
    let lc = target.lower_col;
    for i in 1..=target.row_number() {
        for j in 1..=target.col_number() {
            target.set(lr + i - 1, lc + j - 1, source.get(i as i32, j as i32));
        }
    }
}

/// OCCT AppCont_ContMatrices::MMatrix(classe, M) (the .pxx L67-87) — the
/// Bernstein basis mass matrix of size classe x classe.
pub(crate) fn m_matrix(classe: i32, m: &mut Matrix) {
    assert!(
        classe <= THE_MAX_MMATRIX_CLASSE,
        "Standard_DimensionError: MMatrix: classe > 24"
    );
    let mass = bernstein_mass_matrix(classe);
    fill_with_bounds(m, &mass);
}

/// OCCT AppCont_ContMatrices::InvMMatrix(classe, InvM) (the .pxx L89-106) —
/// the inverse of the Bernstein basis mass matrix.
pub(crate) fn inv_m_matrix(classe: i32, inv_m: &mut Matrix) {
    assert!(
        classe <= THE_MAX_MMATRIX_CLASSE,
        "Standard_DimensionError: InvMMatrix: classe > 24"
    );
    let inverse = bernstein_mass_matrix(classe).inverse();
    fill_with_bounds(inv_m, &inverse);
}

/// OCCT AppCont_ContMatrices::IBPMatrix(classe, IBPMa) (the .pxx L108-128) —
/// the inverse of (MP^T * MP) for pass-point constraints, where MP is the
/// submatrix of M with columns 2..classe-1; result size (classe-2)^2.
pub(crate) fn ibp_matrix(classe: i32, ibp_ma: &mut Matrix) {
    assert!(
        classe <= THE_MAX_CLASSE,
        "Standard_DimensionError: IBPMatrix: classe > 26"
    );
    let inverse = inverse_of_gram_submatrix(classe, 2, classe - 1);
    fill_with_bounds(ibp_ma, &inverse);
}

/// OCCT AppCont_ContMatrices::IBTMatrix(classe, IBTMa) (the .pxx L130-149) —
/// the inverse of (MT^T * MT) for tangency constraints, where MT is the
/// submatrix of M with columns 3..classe-2; result size (classe-4)^2.
pub(crate) fn ibt_matrix(classe: i32, ibt_ma: &mut Matrix) {
    assert!(
        classe <= THE_MAX_CLASSE,
        "Standard_DimensionError: IBTMatrix: classe > 26"
    );
    let inverse = inverse_of_gram_submatrix(classe, 3, classe - 2);
    fill_with_bounds(ibt_ma, &inverse);
}

/// The shared (M_sub^T * M_sub)^-1 computation — the runtime fallback form
/// of AppCont_LeastSquare.cxx L488-499 (math_Matrix MP / IBP1 = MP.Transposed()
/// * MP; IBP = IBP1.Inverse()), applied to the columns [j1..j2] of the mass
/// matrix.
fn inverse_of_gram_submatrix(classe: i32, j1: i32, j2: i32) -> Matrix {
    let mass = bernstein_mass_matrix(classe);
    let nb_col = j2 - j1 + 1;
    let mut mp = Matrix::new(1, classe, 1, nb_col);
    for i in 1..=classe {
        for j in 1..=nb_col {
            mp.set(i, j, mass.get(i, j1 + j - 1));
        }
    }
    mp.transposed().multiplied(&mp).inverse()
}

/// OCCT AppCont_ContMatrices::VBernstein(classe, nbpoints, M) (the .pxx
/// L152-177) — the Bernstein basis values at the Gauss-Legendre points,
/// M(i,j) = B_i(u_j), u_j = (1 - GaussPoint(j)) / 2, filled row-major like
/// the OCCT THE_VBMATRIX_DATA table.
pub(crate) fn v_bernstein(classe: i32, nbpoints: i32, m: &mut Matrix) {
    assert!(
        classe <= THE_MAX_CLASSE,
        "Standard_DimensionError: VBernstein: classe > 26"
    );
    assert!(
        nbpoints <= THE_MAX_NBPOINTS,
        "Standard_DimensionError: VBernstein: nbpoints > 24"
    );

    // The Gauss abscissae (the same OCCT table math::GaussPoints reads).
    // The table columns carry the parameters in ascending order — the same
    // order the least-square constructor produces for myParam (the
    // .cxx L155-170 half reversal of the Gauss sequence): u_j = (1 - G_j)/2
    // on the positive half, mirrored to 1 - u on the second half.
    let mut gauss_p = VecD::new(nbpoints as usize);
    gauss_points(nbpoints as usize, &mut gauss_p);
    let half = (nbpoints + 1) / 2;

    let n = classe - 1;
    let mut vb = Matrix::new(1, classe, 1, nbpoints);
    for i in 1..=classe {
        for j in 1..=nbpoints {
            let u = if j <= half {
                (1.0 - gauss_p.get(j as usize)) * 0.5
            } else {
                1.0 - (1.0 - gauss_p.get((nbpoints + 1 - j) as usize)) * 0.5
            };
            vb.set(i, j, bernstein_value(n, i - 1, u));
        }
    }
    fill_with_bounds(m, &vb);
}

// ---------------------------------------------------------------------------
// OCCT AppCont_LeastSquare (AppCont_LeastSquare.hxx L27-68).
// ---------------------------------------------------------------------------

/// OCCT struct PeriodicityInfo (the .hxx L27-31).
#[derive(Clone, Copy)]
pub struct PeriodicityInfo {
    pub is_periodic: bool,
    pub my_period: f64,
}

/// OCCT AppCont_LeastSquare (the .hxx L33-68) — least square approximation
/// of an AppCont_Function by one Bezier MultiCurve.
pub struct AppContLeastSquare {
    my_scu: MultiCurve,          // OCCT: mySCU
    my_points: Matrix,           // OCCT: myPoints
    my_poles: Matrix,            // OCCT: myPoles
    my_param: RVector,           // OCCT: myParam
    my_vb: Matrix,               // OCCT: myVB
    my_per_info: Vec<PeriodicityInfo>, // OCCT: myPerInfo (Array1 1-based)
    my_done: bool,               // OCCT: myDone
    my_degre: i32,               // OCCT: myDegre
    my_nbdiscret: i32,           // OCCT: myNbdiscret
    my_nb_p: i32,                // OCCT: myNbP
    my_nb_p2d: i32,              // OCCT: myNbP2d
}

impl AppContLeastSquare {
    /// OCCT AppCont_LeastSquare::FixSingleBorderPoint (the .cxx L30-86).
    fn fix_single_border_point(
        &mut self,
        the_ssp: &dyn AppContFunction,
        the_u: f64,
        the_u0: f64,
        the_u1: f64,
        the_fix2d: &mut Vec<DVec2>,
        the_fix: &mut Vec<DVec3>,
    ) {
        let a_max_iter = 15;
        let mut a_tab_p = vec![DVec3::ZERO; (self.my_nb_p.max(1)) as usize];
        let mut a_prev_p = vec![DVec3::ZERO; (self.my_nb_p.max(1)) as usize];
        let mut a_tab_p2d = vec![DVec2::ZERO; (self.my_nb_p2d.max(1)) as usize];
        let mut a_prev_p2d = vec![DVec2::ZERO; (self.my_nb_p2d.max(1)) as usize];
        let a_mult = if (the_u - the_u0) > (the_u1 - the_u) {
            1.0
        } else {
            -1.0
        };
        let a_start_param = the_u;
        let mut a_prev_dist = 1.0f64;
        let mut a_curr_dist = 1.0f64;

        let du = -(the_u1 - the_u0) / 2.0 * a_mult;
        let eps = EPSILON_1;
        let mut dd = du;
        let dec = 0.1;
        for an_iter in 1..a_max_iter {
            dd *= dec;
            let a_curr_param = a_start_param + dd;
            the_ssp.value(a_curr_param, &mut a_tab_p2d, &mut a_tab_p);

            // from second iteration
            if an_iter > 1 {
                a_curr_dist = 0.0;

                let mut i2 = 1;
                for j in 0..self.my_nb_p as usize {
                    a_curr_dist += a_tab_p[j].distance(a_prev_p[j]);
                    i2 += 3;
                }
                for j in 0..self.my_nb_p2d as usize {
                    a_curr_dist += a_tab_p2d[j].distance(a_prev_p2d[j]);
                    i2 += 2;
                }
                let _ = i2; // unused but set for debug

                // from the third iteration
                if an_iter > 2 && a_curr_dist / a_prev_dist > 10.0 {
                    break;
                }
            }
            a_prev_p.copy_from_slice(&a_tab_p);
            a_prev_p2d.copy_from_slice(&a_tab_p2d);
            a_prev_dist = a_curr_dist;
            if a_prev_dist <= eps {
                break;
            }
        }
        *the_fix2d = a_prev_p2d;
        *the_fix = a_prev_p;
    }

    /// OCCT AppCont_LeastSquare::AppCont_LeastSquare(SSP, U0, U1, FirstCons,
    /// LastCons, Deg, myNbPoints) (the .cxx L90-515).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        ssp: &dyn AppContFunction,
        u0: f64,
        u1: f64,
        first_cons: AppParConstraint,
        last_cons: AppParConstraint,
        deg: i32,
        my_nb_points: i32,
    ) -> Self {
        let nb3d = ssp.get_nb_of_3d_points();
        let nb2d = ssp.get_nb_of_2d_points();
        let nbcol = 3 * nb3d + 2 * nb2d;

        let mut ls = AppContLeastSquare {
            my_scu: MultiCurve::new((deg + 1) as usize, 0, 0), // OCCT: mySCU(Deg + 1)
            my_points: Matrix::new(1, my_nb_points, 1, nbcol),
            my_poles: Matrix::new_init(1, deg + 1, 1, nbcol, 0.0),
            my_param: RVector::new(1, my_nb_points),
            my_vb: Matrix::new(1, deg + 1, 1, my_nb_points),
            my_per_info: vec![
                PeriodicityInfo {
                    is_periodic: false,
                    my_period: 0.0,
                };
                (3 * nb3d + 2 * nb2d) as usize
            ],
            my_done: false,
            my_degre: deg,
            my_nbdiscret: my_nb_points,
            my_nb_p: nb3d,
            my_nb_p2d: nb2d,
        };

        // OCCT L104-115.
        let mut my_first_c = first_cons;
        let mut my_last_c = last_cons;
        let classe = deg + 1;
        let cl1 = deg;
        let first_p = 1;
        let last_p = my_nb_points;
        let mut b = Matrix::new_init(1, classe, 1, nbcol, 0.0);
        let mut bdeb = 1;
        let mut bfin = classe;
        let (my_nb_p, my_nb_p2d) = (ls.my_nb_p, ls.my_nb_p2d);

        // OCCT L118-123: the discretisation arrays.
        let mut a_tab_p = vec![DVec3::ZERO; my_nb_p.max(1) as usize];
        let mut a_tab_p2d = vec![DVec2::ZERO; my_nb_p2d.max(1) as usize];
        let mut a_tab_v = vec![DVec3::ZERO; my_nb_p.max(1) as usize];
        let mut a_tab_v2d = vec![DVec2::ZERO; my_nb_p2d.max(1) as usize];

        // OCCT L125-128.
        for a_dim_idx in 1..=(my_nb_p * 3 + my_nb_p2d * 2) {
            let (is_periodic, period) = ssp.period_information(a_dim_idx);
            ls.my_per_info[(a_dim_idx - 1) as usize] = PeriodicityInfo {
                is_periodic,
                my_period: period,
            };
        }

        // OCCT L130-147: the TangencyPoint constraints degrade to PassPoint
        // when the function cannot provide the derivative.
        let mut ok;
        if my_first_c == AppParConstraint::TangencyPoint {
            ok = ssp.d1(u0, &mut a_tab_v2d, &mut a_tab_v);
            if !ok {
                my_first_c = AppParConstraint::PassPoint;
            }
        }

        if my_last_c == AppParConstraint::TangencyPoint {
            ok = ssp.d1(u1, &mut a_tab_v2d, &mut a_tab_v);
            if !ok {
                my_last_c = AppParConstraint::PassPoint;
            }
        }

        // OCCT L149-153: compute control points params on which the
        // approximation will be built (the Gauss integration).
        let mut gauss_p = VecD::new(my_nb_points as usize);
        let mut gauss_w = VecD::new(my_nb_points as usize);
        gauss_points(my_nb_points as usize, &mut gauss_p);
        gauss_weights(my_nb_points as usize, &mut gauss_w);
        let mut the_weights = RVector::new(1, my_nb_points);
        let mut vb_param = RVector::new(1, my_nb_points);
        let du = 0.5 * (u1 - u0);
        for i in first_p..=last_p {
            let u = 0.5 * (u1 + u0) + du * gauss_p.get(i as usize);
            if i <= (my_nb_points + 1) / 2 {
                ls.my_param.set(last_p - i + 1, u);
                vb_param.set(last_p - i + 1, 0.5 * (1.0 + gauss_p.get(i as usize)));
                the_weights.set(last_p - i + 1, 0.5 * gauss_w.get(i as usize));
            } else {
                vb_param.set(i - (my_nb_points + 1) / 2, 0.5 * (1.0 + gauss_p.get(i as usize)));
                ls.my_param.set(i - (my_nb_points + 1) / 2, u);
                the_weights.set(i - (my_nb_points + 1) / 2, 0.5 * gauss_w.get(i as usize));
            }
        }
        let _ = vb_param; // OCCT fills VBParam but never reads it below.

        // OCCT L172-189: compute the control points (the sampled values).
        for i in first_p..=last_p {
            let u = ls.my_param.get(i);
            ssp.value(u, &mut a_tab_p2d, &mut a_tab_p);

            let mut i2 = 1;
            for j in 0..my_nb_p as usize {
                let p = a_tab_p[j];
                ls.my_points.set(i, i2, p.x);
                ls.my_points.set(i, i2 + 1, p.y);
                ls.my_points.set(i, i2 + 2, p.z);
                i2 += 3;
            }
            for j in 0..my_nb_p2d as usize {
                let p = a_tab_p2d[j];
                ls.my_points.set(i, i2, p.x);
                ls.my_points.set(i, i2 + 1, p.y);
                i2 += 2;
            }
        }

        // OCCT L191-221: fix the possible "period jump".
        let a_max_dim = 3 * my_nb_p + 2 * my_nb_p2d;
        for a_dim_idx in 1..=a_max_dim {
            let info = ls.my_per_info[(a_dim_idx - 1) as usize];
            if info.is_periodic
                && (ls.my_points.get(1, a_dim_idx) - ls.my_points.get(2, a_dim_idx)).abs()
                    > info.my_period / 2.01
                && (ls.my_points.get(2, a_dim_idx) - ls.my_points.get(3, a_dim_idx)).abs()
                    < info.my_period / 2.01
            {
                let a_period_mult = if ls.my_points.get(1, a_dim_idx) < ls.my_points.get(2, a_dim_idx)
                {
                    1.0
                } else {
                    -1.0
                };
                let a_new_param =
                    ls.my_points.get(1, a_dim_idx) + a_period_mult * info.my_period;
                ls.my_points.set(1, a_dim_idx, a_new_param);
            }
        }
        for a_pnt_idx in 1..my_nb_points {
            for a_dim_idx in 1..=a_max_dim {
                let info = ls.my_per_info[(a_dim_idx - 1) as usize];
                if info.is_periodic
                    && (ls.my_points.get(a_pnt_idx, a_dim_idx)
                        - ls.my_points.get(a_pnt_idx + 1, a_dim_idx))
                        .abs()
                        > info.my_period / 2.01
                {
                    let a_period_mult = if ls.my_points.get(a_pnt_idx, a_dim_idx)
                        > ls.my_points.get(a_pnt_idx + 1, a_dim_idx)
                    {
                        1.0
                    } else {
                        -1.0
                    };
                    let a_new_param = ls.my_points.get(a_pnt_idx + 1, a_dim_idx)
                        + a_period_mult * info.my_period;
                    ls.my_points.set(a_pnt_idx + 1, a_dim_idx, a_new_param);
                }
            }
        }

        // OCCT L223: AppCont_ContMatrices::VBernstein(classe, myNbPoints, myVB).
        v_bernstein(classe, my_nb_points, &mut ls.my_vb);

        // OCCT L225-243: processing of the right-hand side.
        // (OCCT: NCollection_Array1 tmppoints(1, nbcol) + tmppoints.Init(0.0).)
        let mut tmppoints = RVector::new_init(1, nbcol, 0.0);
        for c in 1..=classe {
            for j in 1..=nbcol {
                tmppoints.set(j, 0.0);
            }
            for i in 1..=my_nb_points {
                let coeff = the_weights.get(i) * ls.my_vb.get(c, i);
                for j in 1..=nbcol {
                    let v = tmppoints.get(j) + ls.my_points.get(i, j) * coeff;
                    tmppoints.set(j, v);
                }
            }
            for k in 1..=nbcol {
                let v = b.get(c, k) + tmppoints.get(k);
                b.set(c, k, v);
            }
        }

        if my_first_c == AppParConstraint::NoConstraint
            && my_last_c == AppParConstraint::NoConstraint
        {
            // OCCT L245-263: direct computation of the poles.
            let mut inv_m = Matrix::new(1, classe, 1, classe);
            inv_m_matrix(classe, &mut inv_m);

            for i in 1..=classe {
                for j in 1..=classe {
                    let i_bij = inv_m.get(i, j);
                    for k in 1..=nbcol {
                        let v = ls.my_poles.get(i, k) + i_bij * b.get(j, k);
                        ls.my_poles.set(i, k, v);
                    }
                }
            }
        } else {
            // OCCT L265-514.
            let mut m = Matrix::new(1, classe, 1, classe);
            m_matrix(classe, &mut m);
            let mut a_fix_p2d = vec![DVec2::ZERO; my_nb_p2d.max(1) as usize];
            let mut a_fix_p = vec![DVec3::ZERO; my_nb_p.max(1) as usize];

            if my_first_c == AppParConstraint::PassPoint
                || my_first_c == AppParConstraint::TangencyPoint
            {
                // OCCT L272-301.
                ssp.value(u0, &mut a_tab_p2d, &mut a_tab_p);
                ls.fix_single_border_point(ssp, u0, u0, u1, &mut a_fix_p2d, &mut a_fix_p);

                let mut i2 = 1;
                for k in 0..my_nb_p as usize {
                    let src = if a_fix_p[k].distance(a_tab_p[k]) > 0.1 {
                        a_fix_p[k]
                    } else {
                        a_tab_p[k]
                    };
                    ls.my_poles.set(1, i2, src.x);
                    ls.my_poles.set(1, i2 + 1, src.y);
                    ls.my_poles.set(1, i2 + 2, src.z);
                    i2 += 3;
                }
                for k in 0..my_nb_p2d as usize {
                    let src = if a_fix_p2d[k].distance(a_tab_p2d[k]) > 0.1 {
                        a_fix_p2d[k]
                    } else {
                        a_tab_p2d[k]
                    };
                    ls.my_poles.set(1, i2, src.x);
                    ls.my_poles.set(1, i2 + 1, src.y);
                    i2 += 2;
                }

                // OCCT L303-312.
                for a_dim_idx in 1..=a_max_dim {
                    let info = ls.my_per_info[(a_dim_idx - 1) as usize];
                    if info.is_periodic
                        && (ls.my_poles.get(1, a_dim_idx) - ls.my_points.get(1, a_dim_idx)).abs()
                            > info.my_period / 2.01
                    {
                        let a_mult = if ls.my_poles.get(1, a_dim_idx)
                            < ls.my_points.get(1, a_dim_idx)
                        {
                            1.0
                        } else {
                            -1.0
                        };
                        let v = ls.my_poles.get(1, a_dim_idx) + a_mult * info.my_period;
                        ls.my_poles.set(1, a_dim_idx, v);
                    }
                }
            }

            if my_last_c == AppParConstraint::PassPoint
                || my_last_c == AppParConstraint::TangencyPoint
            {
                // OCCT L315-356.
                ssp.value(u1, &mut a_tab_p2d, &mut a_tab_p);
                ls.fix_single_border_point(ssp, u1, u0, u1, &mut a_fix_p2d, &mut a_fix_p);

                let mut i2 = 1;
                for k in 0..my_nb_p as usize {
                    let src = if a_fix_p[k].distance(a_tab_p[k]) > 0.1 {
                        a_fix_p[k]
                    } else {
                        a_tab_p[k]
                    };
                    ls.my_poles.set(classe, i2, src.x);
                    ls.my_poles.set(classe, i2 + 1, src.y);
                    ls.my_poles.set(classe, i2 + 2, src.z);
                    i2 += 3;
                }
                for k in 0..my_nb_p2d as usize {
                    let src = if a_fix_p2d[k].distance(a_tab_p2d[k]) > 0.1 {
                        a_fix_p2d[k]
                    } else {
                        a_tab_p2d[k]
                    };
                    ls.my_poles.set(classe, i2, src.x);
                    ls.my_poles.set(classe, i2 + 1, src.y);
                    i2 += 2;
                }

                // OCCT L346-355 — the OCCT loop bound is the literal 2 (not
                // aMaxDim); the quirk is kept as-is.
                for a_dim_idx in 1..=2 {
                    let info = ls.my_per_info[(a_dim_idx - 1) as usize];
                    if info.is_periodic
                        && (ls.my_poles.get(classe, a_dim_idx)
                            - ls.my_points.get(my_nb_points, a_dim_idx))
                            .abs()
                            > info.my_period / 2.01
                    {
                        let a_mult = if ls.my_poles.get(classe, a_dim_idx)
                            < ls.my_points.get(my_nb_points, a_dim_idx)
                        {
                            1.0
                        } else {
                            -1.0
                        };
                        let v =
                            ls.my_poles.get(classe, a_dim_idx) + a_mult * info.my_period;
                        ls.my_poles.set(classe, a_dim_idx, v);
                    }
                }
            }

            if my_first_c == AppParConstraint::PassPoint {
                // OCCT L358-370.
                bdeb = 2;
                for i in 1..=classe {
                    let coeff = m.get(i, 1);
                    for k in 1..=nbcol {
                        let v = b.get(i, k) - ls.my_poles.get(1, k) * coeff;
                        b.set(i, k, v);
                    }
                }
            }

            if my_last_c == AppParConstraint::PassPoint {
                // OCCT L372-383.
                bfin = cl1;
                for i in 1..=classe {
                    let coeff = m.get(i, classe);
                    for k in 1..=nbcol {
                        let v = b.get(i, k) - ls.my_poles.get(classe, k) * coeff;
                        b.set(i, k, v);
                    }
                }
            }

            if my_first_c == AppParConstraint::TangencyPoint {
                // OCCT L385-419: fix the second pole.
                bdeb = 3;
                ssp.d1(u0, &mut a_tab_v2d, &mut a_tab_v);

                let mut i2 = 1;
                let coeff = (u1 - u0) / ls.my_degre as f64;
                for k in 0..my_nb_p as usize {
                    let i2plus1 = i2 + 1;
                    let i2plus2 = i2 + 2;
                    ls.my_poles.set(2, i2, ls.my_poles.get(1, i2) + a_tab_v[k].x * coeff);
                    ls.my_poles.set(
                        2,
                        i2plus1,
                        ls.my_poles.get(1, i2plus1) + a_tab_v[k].y * coeff,
                    );
                    ls.my_poles.set(
                        2,
                        i2plus2,
                        ls.my_poles.get(1, i2plus2) + a_tab_v[k].z * coeff,
                    );
                    i2 += 3;
                }
                for k in 0..my_nb_p2d as usize {
                    let i2plus1 = i2 + 1;
                    ls.my_poles
                        .set(2, i2, ls.my_poles.get(1, i2) + a_tab_v2d[k].x * coeff);
                    ls.my_poles.set(
                        2,
                        i2plus1,
                        ls.my_poles.get(1, i2plus1) + a_tab_v2d[k].y * coeff,
                    );
                    i2 += 2;
                }

                for i in 1..=classe {
                    let coeff = m.get(i, 1);
                    let coeff2 = m.get(i, 2);
                    for k in 1..=nbcol {
                        let v = b.get(i, k)
                            - (ls.my_poles.get(1, k) * coeff + ls.my_poles.get(2, k) * coeff2);
                        b.set(i, k, v);
                    }
                }
            }

            if my_last_c == AppParConstraint::TangencyPoint {
                // OCCT L421-453.
                bfin = classe - 2;
                ssp.d1(u1, &mut a_tab_v2d, &mut a_tab_v);
                let mut i2 = 1;
                let coeff = (u1 - u0) / ls.my_degre as f64;
                for k in 0..my_nb_p as usize {
                    let i2plus1 = i2 + 1;
                    let i2plus2 = i2 + 2;
                    ls.my_poles.set(
                        cl1,
                        i2,
                        ls.my_poles.get(classe, i2) - a_tab_v[k].x * coeff,
                    );
                    ls.my_poles.set(
                        cl1,
                        i2plus1,
                        ls.my_poles.get(classe, i2plus1) - a_tab_v[k].y * coeff,
                    );
                    ls.my_poles.set(
                        cl1,
                        i2plus2,
                        ls.my_poles.get(classe, i2plus2) - a_tab_v[k].z * coeff,
                    );
                    i2 += 3;
                }
                for k in 0..my_nb_p2d as usize {
                    let i2plus1 = i2 + 1;
                    ls.my_poles.set(
                        cl1,
                        i2,
                        ls.my_poles.get(classe, i2) - a_tab_v2d[k].x * coeff,
                    );
                    ls.my_poles.set(
                        cl1,
                        i2plus1,
                        ls.my_poles.get(classe, i2plus1) - a_tab_v2d[k].y * coeff,
                    );
                    i2 += 2;
                }

                for i in 1..=classe {
                    let coeff = m.get(i, classe);
                    let coeff2 = m.get(i, cl1);
                    for k in 1..=nbcol {
                        let v = b.get(i, k)
                            - (ls.my_poles.get(classe, k) * coeff
                                + ls.my_poles.get(cl1, k) * coeff2);
                        b.set(i, k, v);
                    }
                }
            }

            if bdeb <= bfin {
                // OCCT L455-513.
                let mut b2 = Matrix::new_init(bdeb, bfin, 1, b.upper_col(), 0.0);

                for i in bdeb..=bfin {
                    for j in 1..=classe {
                        let coeff = m.get(i, j);
                        for k in 1..=nbcol {
                            let v = b2.get(i, k) + b.get(j, k) * coeff;
                            b2.set(i, k, v);
                        }
                    }
                }

                // OCCT L471-499: resolution.
                let ibp;

                // In IBPMatrix and IBTMatrix only results for a degree
                // less than or equal to 26 are stored (for now at least).
                if bdeb == 2 && bfin == classe - 1 && classe <= 26 {
                    let mut m_ibp = Matrix::new(bdeb, bfin, bdeb, bfin);
                    ibp_matrix(classe, &mut m_ibp);
                    ibp = m_ibp;
                } else if bdeb == 3 && bfin == classe - 2 && classe <= 26 {
                    let mut m_ibt = Matrix::new(bdeb, bfin, bdeb, bfin);
                    ibt_matrix(classe, &mut m_ibt);
                    ibp = m_ibt;
                } else {
                    let mut mp = Matrix::new(1, classe, bdeb, bfin);
                    for i in 1..=classe {
                        for j in bdeb..=bfin {
                            mp.set(i, j, m.get(i, j));
                        }
                    }
                    // OCCT: IBP1 = MP.Transposed() * MP; IBP = IBP1.Inverse()
                    // — the assigned matrix keeps the (bdeb..bfin) bounds.
                    let inverted = mp.transposed().multiplied(&mp).inverse();
                    let mut m_fallback = Matrix::new(bdeb, bfin, bdeb, bfin);
                    fill_with_bounds(&mut m_fallback, &inverted);
                    ibp = m_fallback;
                }

                ls.my_done = true;
                for i in bdeb..=bfin {
                    for j in bdeb..=bfin {
                        let i_bpij = ibp.get(i, j);
                        for k in 1..=nbcol {
                            let v = ls.my_poles.get(i, k) + i_bpij * b2.get(j, k);
                            ls.my_poles.set(i, k, v);
                        }
                    }
                }
            }
        }

        ls
    }

    /// OCCT AppCont_LeastSquare::Value() (the .cxx L519-547).
    pub fn value(&mut self) -> MultiCurve {
        let ideb = 1;
        let ifin = self.my_degre + 1;
        let my_nb_p = self.my_nb_p as usize;
        let my_nb_p2d = self.my_nb_p2d as usize;

        // Store the result in the corresponding curves.
        for i in ideb..=ifin {
            let mut j2 = 1usize;
            let mut m_pole = MultiPoint::new(my_nb_p, my_nb_p2d);
            for j in 1..=my_nb_p {
                let p = DVec3::new(
                    self.my_poles.get(i, j2 as i32),
                    self.my_poles.get(i, (j2 + 1) as i32),
                    self.my_poles.get(i, (j2 + 2) as i32),
                );
                m_pole.set_point(j, p);
                j2 += 3;
            }
            for j in (my_nb_p + 1)..=(my_nb_p + my_nb_p2d) {
                let p = DVec2::new(
                    self.my_poles.get(i, j2 as i32),
                    self.my_poles.get(i, (j2 + 1) as i32),
                );
                m_pole.set_point2d(j, p);
                j2 += 2;
            }
            self.my_scu.set_value(i as usize, m_pole);
        }
        self.my_scu.clone()
    }

    /// OCCT AppCont_LeastSquare::Error(F, MaxE3d, MaxE2d) (the .cxx L551-607).
    pub fn error(&self) -> (f64, f64, f64) {
        let classe = self.my_degre + 1;
        let mut err3d = 0.0f64;
        let mut err2d = 0.0f64;
        let ncol = self.my_points.col_number();

        let mut my_points = Matrix::new(1, self.my_nbdiscret, 1, ncol);
        for i in 1..=self.my_nbdiscret {
            for j in 1..=ncol {
                my_points.set(i, j, self.my_points.get(i, j));
            }
        }

        let mut max_e3d = 0.0f64;
        let mut max_e2d = 0.0f64;
        let mut f = 0.0f64;

        let mut tmppoles = RVector::new(1, ncol);

        for c in 1..=classe {
            for k in 1..=ncol {
                tmppoles.set(k, self.my_poles.get(c, k));
            }
            for i in 1..=self.my_nbdiscret {
                let coeff = self.my_vb.get(c, i);
                for j in 1..=ncol {
                    let v = my_points.get(i, j) - tmppoles.get(j) * coeff;
                    my_points.set(i, j, v);
                }
            }
        }

        for i in 1..=self.my_nbdiscret {
            let mut i2 = 1;
            for _j in 1..=self.my_nb_p {
                let e1 = my_points.get(i, i2);
                let e2 = my_points.get(i, i2 + 1);
                let e3 = my_points.get(i, i2 + 2);
                err3d = e1 * e1 + e2 * e2 + e3 * e3;
                max_e3d = max_e3d.max(err3d);
                f += err3d;
                i2 += 3;
            }
            for _j in 1..=self.my_nb_p2d {
                let e1 = my_points.get(i, i2);
                let e2 = my_points.get(i, i2 + 1);
                err2d = e1 * e1 + e2 * e2;
                max_e2d = max_e2d.max(err2d);
                f += err2d;
                i2 += 2;
            }
        }

        max_e3d = max_e3d.sqrt();
        max_e2d = max_e2d.sqrt();
        (f, max_e3d, max_e2d)
    }

    /// OCCT AppCont_LeastSquare::IsDone() (the .cxx L611-614).
    pub fn is_done(&self) -> bool {
        self.my_done
    }
}

// The Precision::Confusion read of the Approx engine (kept here so the
// app_cont unit tests exercise the same constants the least square uses).
const _: f64 = CONFUSION;

/// OCCT Standard_Real RealLast() re-host for the module (Standard_Real.hxx
/// L93-96) — used by the Approx_ComputeCLine engine.
pub(crate) fn real_last() -> f64 {
    REAL_LAST
}

/// OCCT standard_epsilon on an arbitrary value — re-exported for the
/// Approx_FitAndDivide module tests.
#[allow(dead_code)]
pub(crate) fn epsilon_of(v: f64) -> f64 {
    standard_epsilon(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    // The expected values below are literal entries of the OCCT generated
    // tables (AppCont_ContMatrices_{VB,BB,InvM,IBP,IBT}.pxx), extracted at
    // the (classe=9, nbpoints=17) block — the degree-8 path of the
    // BiTgte_Blend MakeCurve consumer.

    /// VBernstein(9, 17) vs the OCCT THE_VBMATRIX_DATA block at offset
    /// 300*((9-1)*9/2 - 1) + 9*16*17/2 = 11724.
    #[test]
    fn v_bernstein_matches_occt_table_classe9_nbpoints17() {
        let mut m = Matrix::new(1, 9, 1, 17);
        v_bernstein(9, 17, &mut m);
        // OCCT VB[11724..11729] — the first row B_0(u_j), j = 1..5.
        let expected_first_row = [
            0.96291782758848,
            0.818917795111741,
            0.61018955814276,
            0.396345524976371,
            0.222713778885278,
        ];
        for (j, want) in expected_first_row.iter().enumerate() {
            let got = m.get(1, j as i32 + 1);
            assert!(
                (got - want).abs() < 1e-12,
                "VB(1,{}) = {} vs OCCT {}",
                j + 1,
                got,
                want
            );
        }
        // OCCT VB[11724 + 152] — the last entry B_8(u_17) (symmetry with the
        // first entry of the mirrored Gauss point).
        let got = m.get(9, 17);
        assert!((got - 0.96291782758848).abs() < 1e-12, "VB(9,17) = {}", got);
    }

    /// MMatrix(9) vs the OCCT THE_BBMATRIX_DATA block at offset
    /// 4+9+16+25+36+49+64 = 203 (1/9, 1/18, ... exact rationals).
    #[test]
    fn m_matrix_matches_occt_table_classe9() {
        let mut m = Matrix::new(1, 9, 1, 9);
        m_matrix(9, &mut m);
        let expected = [
            0.05882352941176470588235294, // BB[203] = M(1,1) = 1/9
            0.02941176470588235294117647, // BB[204] = M(1,2)
            0.01372549019607843137254902, // BB[205] = M(1,3)
            0.005882352941176470588235294, // BB[206] = M(1,4)
            0.002262443438914027149321267, // BB[207] = M(1,5)
            0.000754147812971342383107089, // BB[208] = M(1,6)
        ];
        for (j, want) in expected.iter().enumerate() {
            let got = m.get(1, j as i32 + 1);
            assert!(
                (got - want).abs() < 1e-15,
                "M(1,{}) = {} vs OCCT {}",
                j + 1,
                got,
                want
            );
        }
        // OCCT BB[203 + 80] = M(9,9) = 1/9.
        let got = m.get(9, 9);
        assert!((got - 0.05882352941176470588235294).abs() < 1e-15);
    }

    /// InvMMatrix(9) vs the OCCT THE_INV_MMATRIX_DATA block at offset 203 —
    /// the first row is exactly integral: 81, -324, 756, -1134, 1134, -756.
    #[test]
    fn inv_m_matrix_matches_occt_table_classe9() {
        let mut m = Matrix::new(1, 9, 1, 9);
        inv_m_matrix(9, &mut m);
        let expected = [81.0, -324.0, 756.0, -1134.0, 1134.0, -756.0];
        for (j, want) in expected.iter().enumerate() {
            let got = m.get(1, j as i32 + 1);
            assert!(
                (got - want).abs() < 1e-6,
                "InvM(1,{}) = {} vs OCCT {}",
                j + 1,
                got,
                want
            );
        }
    }

    /// IBPMatrix(9) vs the OCCT THE_IBPMATRIX_DATA block at offset 91 —
    /// the (classe-2)^2 = 7x7 inverse for the PassPoint path.
    #[test]
    fn ibp_matrix_matches_occt_table_classe9() {
        let mut m = Matrix::new(1, 7, 1, 7);
        ibp_matrix(9, &mut m);
        let expected = [
            379_424.357_666_301_989_822_182_4,
            -2.100_567_756_809_113_573_889_156e6,
            4.937_601_766_796_749_036_011_036e6,
            -6.271_766_132_734_029_342_492_488e6,
            4.529_063_696_431_671_303_622_051e6,
            -1.761_378_397_072_727_762_621_633e6,
        ];
        for (j, want) in expected.iter().enumerate() {
            let got = m.get(1, j as i32 + 1);
            assert!(
                (got - want).abs() <= 1e-6 * want.abs(),
                "IBP(1,{}) = {} vs OCCT {}",
                j + 1,
                got,
                want
            );
        }
    }

    /// IBTMatrix(9) vs the OCCT THE_IBTMATRIX_DATA block at offset 30 —
    /// the (classe-4)^2 = 5x5 inverse for the TangencyPoint path.
    #[test]
    fn ibt_matrix_matches_occt_table_classe9() {
        let mut m = Matrix::new(1, 5, 1, 5);
        ibt_matrix(9, &mut m);
        let expected = [
            131_995.562_952_392_883_982_868_6,
            -459_851.698_594_009_545_424_292,
            633_555.123_965_057_551_337_059,
            -405_583.140_087_099_050_069_762_9,
            101_913.355_235_452_463_620_874_2,
        ];
        for (j, want) in expected.iter().enumerate() {
            let got = m.get(1, j as i32 + 1);
            assert!(
                (got - want).abs() <= 1e-6 * want.abs(),
                "IBT(1,{}) = {} vs OCCT {}",
                j + 1,
                got,
                want
            );
        }
    }

    /// OCCT Standard::Epsilon sanity (Standard_Real.hxx L242-250): the
    /// neighbour distance of 1.0 is DBL_EPSILON.
    #[test]
    fn standard_epsilon_matches_occt() {
        assert!((epsilon_of(1.0) - f64::EPSILON).abs() < 1e-30);
        assert!(epsilon_of(16.0) > 0.0);
    }
}
