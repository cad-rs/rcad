//! OCCT BSplCLib (TKMath) — knot sequences, interpolation, evaluation,
//! `MovePointAndTangent` and degree elevation.
//!
//! 1:1 translation of the pieces of `BSplCLib.cxx` / `BSplCLib_2.cxx`
//! consumed by the helix pipeline (`Convert_CompPolynomialToPoles`,
//! `Geom_BSplineCurve::MovePointAndTangent`,
//! `Geom_BSplineCurve::IncreaseDegree`): Hunt, FlatIndex, NbPoles,
//! KnotSequenceLength, KnotSequence, LocateParameter (used overloads),
//! BuildSchoenbergPoints, EvalBsplineBasis, BuildBSpMatrix,
//! FactorBandedMatrix, SolveBandedSystem, Interpolate, Eval (non-rational and
//! homogeneous flat-knots overloads), MovePointAndTangent, BuildKnots,
//! BuildBoor, BoorIndex, GetPole, BoorScheme, Copy, InsertKnots,
//! IncreaseDegreeCountKnots, IncreaseDegree, plus
//! `PLib::RationalDerivatives` (PLib.cxx).
//!
//! Index convention: OCCT arrays are 1-based.  Helpers [`at`] / [`ati`] /
//! [`set_at`] keep the OCCT index arithmetic verbatim (subtracting one at the
//! slice access).  `BandMatrix` stands in for the band-shaped `math_Matrix`.

use super::plib::eval_polynomial_flat;

/// OCCT BSplCLib::FirstUKnotIndex(Degree, Mults) over (knots, mults) arrays.
pub fn first_uknot_index_mults(degree: usize, mults: &[i32]) -> i32 {
    let mut index = 1i32;
    let mut sigma = ati(mults, index);
    while sigma <= degree as i32 {
        index += 1;
        sigma += ati(mults, index);
    }
    index
}

/// OCCT BSplCLib::LastUKnotIndex(Degree, Mults) over (knots, mults) arrays.
pub fn last_uknot_index_mults(degree: usize, mults: &[i32]) -> i32 {
    let mut index = mults.len() as i32;
    let mut sigma = ati(mults, index);
    while sigma <= degree as i32 {
        index -= 1;
        sigma += ati(mults, index);
    }
    index
}

/// OCCT gp::Resolution().
pub const GP_RESOLUTION: f64 = 1.0e-12;
/// OCCT Geom_BSplineCurve::MaxDegree() == BSplCLib::MaxDegree().
pub const BSPLIB_MAX_DEGREE: usize = 25;

/// Read `arr[i]` with OCCT 1-based index semantics.
#[inline]
pub fn at(arr: &[f64], i: i32) -> f64 {
    arr[(i - 1) as usize]
}
/// Read an i32 array with OCCT 1-based index semantics.
#[inline]
pub fn ati(arr: &[i32], i: i32) -> i32 {
    arr[(i - 1) as usize]
}
/// Write `arr[i] = v` with OCCT 1-based index semantics.
#[inline]
pub fn set_at(arr: &mut [f64], i: i32, v: f64) {
    arr[(i - 1) as usize] = v;
}
/// Write `arr[i] = v` (i32) with OCCT 1-based index semantics.
#[inline]
pub fn set_ati(arr: &mut [i32], i: i32, v: i32) {
    arr[(i - 1) as usize] = v;
}

/// OCCT Epsilon(1.) — relative machine epsilon.
#[inline]
fn epsilon() -> f64 {
    f64::EPSILON
}

/// OCCT RealSmall().
#[inline]
fn real_small() -> f64 {
    2.2250738585072014e-308 // DBL_MIN
}

/// OCCT BSplCLib::Hunt (BSplCLib.cxx L74-102) — dichotomy; on exit
/// `x_pos` is the largest index with `array[x_pos] < x` (1-based; 0 if x <=
/// array[1], len+1 if x >= array[len]).
pub fn hunt(array: &[f64], x: f64, x_pos: &mut i32) {
    let lower = 1i32;
    let upper = array.len() as i32;
    if array[0] > x {
        *x_pos = lower - 1;
        return;
    } else if array[(upper - 1) as usize] < x {
        *x_pos = upper + 1;
        return;
    }

    *x_pos = lower;
    if upper - lower <= 0 {
        return;
    }

    let mut hi = upper;
    while hi - *x_pos != 1 {
        let mid = (hi + *x_pos) / 2;
        if at(array, mid) < x {
            *x_pos = mid;
        } else {
            hi = mid;
        }
    }
}

/// OCCT BSplCLib::FlatIndex.
pub fn flat_index(degree: usize, index: usize, mults: &[i32], periodic: bool) -> usize {
    let degree = degree as i32;
    let mut index = index as i32;
    for i in 2..=(index) {
        index += ati(mults, i) - 1;
    }
    if periodic {
        index += degree;
    } else {
        index += ati(mults, 1) - 1;
    }
    index as usize
}

/// OCCT BSplCLib::NbPoles(Degree, Periodic, Mults) (BSplCLib.cxx L392-451).
pub fn nb_poles(degree: usize, periodic: bool, mults: &[i32]) -> usize {
    let l = mults.len() as i32;
    let mf = ati(mults, 1);
    let ml = ati(mults, l);
    if mf <= 0 {
        return 0;
    }
    if ml <= 0 {
        return 0;
    }
    let degree = degree as i32;
    let mut sigma;
    if periodic {
        if mf > degree {
            return 0;
        }
        if ml > degree {
            return 0;
        }
        if mf != ml {
            return 0;
        }
        sigma = mf;
    } else {
        let deg1 = degree + 1;
        if mf > deg1 {
            return 0;
        }
        if ml > deg1 {
            return 0;
        }
        sigma = mf + ml - deg1;
    }

    for i in 1..(l - 1) {
        // OCCT: pmu[i] = Mults(Mults.Lower() + i), i = 1..(Upper - Lower - 1).
        let m = ati(mults, 1 + i);
        if m <= 0 {
            return 0;
        }
        if m > degree {
            return 0;
        }
        sigma += m;
    }
    sigma.max(0) as usize
}

/// OCCT BSplCLib::KnotSequenceLength(Mults, Degree, Periodic).
pub fn knot_sequence_length(mults: &[i32], degree: usize, periodic: bool) -> usize {
    let mut l = 0i32;
    let m_upper = mults.len() as i32;
    for i in 1..=m_upper {
        l += ati(mults, i);
    }
    if periodic {
        l += 2 * (degree as i32 + 1 - ati(mults, 1));
    }
    l as usize
}

/// OCCT BSplCLib::KnotSequence(Knots, Mults, Degree, Periodic, KnotSeq)
/// (BSplCLib.cxx L488-547).  `knot_seq` must be sized by
/// [`knot_sequence_length`].
pub fn knot_sequence(
    knots: &[f64],
    mults: &[i32],
    degree: usize,
    periodic: bool,
    knot_seq: &mut [f64],
) {
    let m1 = degree as i32 + 1 - ati(mults, 1); // for periodic
    let k_upper = knots.len() as i32;
    let mut index = if periodic { m1 + 1 } else { 1 };

    for i in 1..=k_upper {
        let mult = ati(mults, i);
        let k = at(knots, i);
        for _ in 1..=mult {
            set_at(knot_seq, index, k);
            index += 1;
        }
    }
    if periodic {
        let period = at(knots, k_upper) - at(knots, 1);
        let mut m = 1;
        let mut j = k_upper - 1;
        for i in (1..=m1).rev() {
            set_at(knot_seq, i, at(knots, j) - period);
            m += 1;
            if m > ati(mults, j) {
                j -= 1;
                m = 1;
            }
        }
        let mut m = 1;
        let mut j = 2i32;
        let upper = knot_seq.len() as i32;
        for i in index..=upper {
            set_at(knot_seq, i, at(knots, j) + period);
            m += 1;
            if m > ati(mults, j) {
                j += 1;
                m = 1;
            }
        }
    }
}

/// OCCT BSplCLib::LocateParameter(Knots, U, IsPeriodic, FromK1, ToK2,
/// KnotIndex, NewU, UFirst, ULast) (BSplCLib.cxx L218-316).
pub fn locate_parameter_main(
    knots: &[f64],
    u: f64,
    is_periodic: bool,
    from_k1: i32,
    to_k2: i32,
    knot_index: &mut i32,
    new_u: &mut f64,
    ufirst: f64,
    ulast: f64,
) {
    let (first, last) = if from_k1 < to_k2 {
        (from_k1, to_k2)
    } else {
        (to_k2, from_k1)
    };
    let last1 = last - 1;
    *new_u = u;
    if is_periodic && (*new_u < ufirst || *new_u > ulast) {
        // ElCLib::InPeriod(NewU, UFirst, ULast).
        let period = ulast - ufirst;
        if period.abs() > f64::EPSILON {
            let mut v = *new_u - ufirst;
            v -= (v / period).floor() * period;
            *new_u = ufirst + v;
        }
    }

    hunt(knots, *new_u, knot_index);

    let k_upper = knots.len() as i32;
    let eps = epsilon() * at(knots, k_upper).abs().min(u.abs());

    if *knot_index < k_upper {
        let mut val = *new_u - at(knots, *knot_index + 1);
        if val < 0.0 {
            val = -val;
        }
        // <= to be coherent with Segment where Eps corresponds to a bit of error.
        if val <= eps {
            *knot_index += 1;
        }
    }
    if *knot_index < first {
        *knot_index = first;
    }
    if *knot_index > last1 {
        *knot_index = last1;
    }

    if *knot_index != last1 {
        let mut k1 = at(knots, *knot_index);
        let mut k2 = at(knots, *knot_index + 1);
        let mut val = k2 - k1;
        if val < 0.0 {
            val = -val;
        }
        while val <= eps {
            *knot_index += 1;
            if *knot_index >= k_upper {
                break;
            }
            k1 = k2;
            k2 = at(knots, *knot_index + 1);
            val = k2 - k1;
            if val < 0.0 {
                val = -val;
            }
        }
    }
}

/// OCCT BSplCLib::LocateParameter(Degree, Knots, U, IsPeriodic, FromK1, ToK2,
/// KnotIndex, NewU) — flat-knots overload (BSplCLib.cxx L189-217).
pub fn locate_parameter_flat(
    degree: usize,
    knots: &[f64],
    u: f64,
    is_periodic: bool,
    from_k1: i32,
    to_k2: i32,
    knot_index: &mut i32,
    new_u: &mut f64,
) {
    let degree = degree as i32;
    let k_upper = knots.len() as i32;
    if is_periodic {
        locate_parameter_main(
            knots,
            u,
            is_periodic,
            from_k1,
            to_k2,
            knot_index,
            new_u,
            at(knots, 1 + degree),
            at(knots, k_upper - degree),
        );
    } else {
        locate_parameter_main(knots, u, is_periodic, from_k1, to_k2, knot_index, new_u, 0.0, 1.0);
    }
}

/// OCCT BSplCLib::BuildSchoenbergPoints (BSplCLib.cxx L3331-3350).
pub fn build_schoenberg_points(degree: usize, flat_knots: &[f64], parameters: &mut [f64]) {
    let inverse = 1.0 / degree as f64;
    let n = parameters.len() as i32;
    for ii in 1..=n {
        let mut sum = 0.0f64;
        for jj in 1..=(degree as i32) {
            sum += at(flat_knots, jj + ii);
        }
        parameters[(ii - 1) as usize] = sum * inverse;
    }
}

/// Band matrix helper: rows `1..=n_rows`, cols `1..=n_cols`, row-major flat
/// storage standing in for the OCCT `math_Matrix`.
pub struct BandMatrix {
    pub n_rows: i32,
    pub n_cols: i32,
    pub data: Vec<f64>,
}

impl BandMatrix {
    pub fn new(n_rows: i32, n_cols: i32) -> Self {
        BandMatrix {
            n_rows,
            n_cols,
            data: vec![0.0; (n_rows * n_cols) as usize],
        }
    }
    #[inline]
    pub fn get(&self, i: i32, j: i32) -> f64 {
        self.data[((i - 1) * self.n_cols + (j - 1)) as usize]
    }
    #[inline]
    pub fn set(&mut self, i: i32, j: i32, v: f64) {
        self.data[((i - 1) * self.n_cols + (j - 1)) as usize] = v;
    }
    pub fn init(&mut self, v: f64) {
        self.data.fill(v);
    }
}

/// OCCT BSplCLib::EvalBsplineBasis (BSplCLib_2.cxx L429-...).
/// `bspline_basis` must have at least `local_request + 1` rows and `order`
/// columns; row `r` receives derivative `r-1` of the non-zero basis.
pub fn eval_bspline_basis(
    derivative_request: i32,
    order: i32,
    flat_knots: &[f64],
    parameter: f64,
    first_non_zero_bspline_index: &mut i32,
    bspline_basis: &mut BandMatrix,
    is_periodic: bool,
) -> i32 {
    let mut local_request = derivative_request;
    if derivative_request >= order {
        local_request = order - 1;
    }

    let k_upper = flat_knots.len() as i32;
    let a_num_poles = k_upper - order; // OCCT: Upper() - Lower() + 1 - Order (Lower == 1).
    let mut ii = 0i32;
    let mut new_param = 0.0f64;
    locate_parameter_flat(
        (order - 1) as usize,
        flat_knots,
        parameter,
        is_periodic,
        order,
        a_num_poles + 1,
        &mut ii,
        &mut new_param,
    );

    *first_non_zero_bspline_index = ii - order + 1;

    let resolution = GP_RESOLUTION;
    ii -= 1; // rebase: `ii` now indexes the knots 0-based (OCCT: ii -= Lower())

    let ncols = bspline_basis.n_cols;

    // Flat access: aBasisData[k-1] (the plain de Boor value array) is the
    // flat index k-1 (row 1).  Derivative rows start at multiples of ncols.
    let basis = &mut bspline_basis.data;

    // --- plain de Boor values up to derivative order's first phase.
    basis[0] = 1.0;

    for qq in 2..=(order - local_request) {
        basis[(qq - 1) as usize] = 0.0;
        for pp in 1..=(qq - 1) {
            let scale = flat_knots[(ii + pp) as usize] - flat_knots[(ii - qq + pp + 1) as usize];
            if scale.abs() < resolution {
                return 2;
            }
            let factor = (parameter - flat_knots[(ii - qq + pp + 1) as usize]) / scale;
            let saved = factor * basis[(pp - 1) as usize];
            basis[(pp - 1) as usize] *= 1.0 - factor;
            basis[(pp - 1) as usize] += basis[(qq - 1) as usize];
            basis[(qq - 1) as usize] = saved;
        }
    }

    for qq in (order - local_request + 1)..=order {
        for pp in 1..=(qq - 1) {
            basis[((order - qq + 1) * ncols + (pp - 1)) as usize] = basis[(pp - 1) as usize];
        }
        basis[(qq - 1) as usize] = 0.0;
        for ss in (order - local_request + 1)..=qq {
            basis[((order - ss + 1) * ncols + (qq - 1)) as usize] = 0.0;
        }
        for pp in 1..=(qq - 1) {
            let scale = flat_knots[(ii + pp) as usize] - flat_knots[(ii - qq + pp + 1) as usize];
            if scale.abs() < resolution {
                return 2;
            }
            let inverse = 1.0 / scale;
            let factor = (parameter - flat_knots[(ii - qq + pp + 1) as usize]) * inverse;
            let mut saved = factor * basis[(pp - 1) as usize];
            basis[(pp - 1) as usize] *= 1.0 - factor;
            basis[(pp - 1) as usize] += basis[(qq - 1) as usize];
            basis[(qq - 1) as usize] = saved;
            let local_inverse = (qq - 1) as f64 * inverse;
            for ss in (order - local_request + 1)..=qq {
                let row_s = ((order - ss + 1) * ncols) as usize;
                saved = local_inverse * basis[row_s + (pp - 1) as usize];
                basis[row_s + (pp - 1) as usize] *= -local_inverse;
                basis[row_s + (pp - 1) as usize] += basis[row_s + (qq - 1) as usize];
                basis[row_s + (qq - 1) as usize] = saved;
            }
        }
    }

    0
}

/// OCCT BSplCLib::BuildBSpMatrix (BSplCLib_2.cxx L327-...).
pub fn build_bsp_matrix(
    parameters: &[f64],
    contact_order_array: &[i32],
    flat_knots: &[f64],
    degree: usize,
    matrix: &mut BandMatrix,
    upper_band_width: &mut i32,
    lower_band_width: &mut i32,
) -> i32 {
    let a_max_order = BSPLIB_MAX_DEGREE as i32 + 1;
    let order = degree as i32 + 1;
    *upper_band_width = degree as i32;
    *lower_band_width = degree as i32;

    if matrix.n_rows != parameters.len() as i32 || matrix.n_cols != 2 * degree as i32 + 1 {
        return 1;
    }

    let mut a_bspline_basis = BandMatrix::new(a_max_order, a_max_order);
    matrix.init(0.0);

    for i in 1..=parameters.len() as i32 {
        let mut first_non_zero_index = 0i32;
        let contact_order = contact_order_array[(i - 1) as usize];
        let an_error_code = eval_bspline_basis(
            contact_order,
            order,
            flat_knots,
            parameters[(i - 1) as usize],
            &mut first_non_zero_index,
            &mut a_bspline_basis,
            false,
        );
        if an_error_code != 0 {
            return 2;
        }

        let an_index = *lower_band_width + 1 + first_non_zero_index - i;
        // aBasisSrc = aBasisBuf + aContactOrder * aMaxOrder → row
        // (contact_order + 1), entries j = 0..order → columns 1..=order.
        for j in 0..order {
            let v = a_bspline_basis.get(contact_order + 1, j + 1);
            matrix.set(i, an_index + j, v);
        }
    }

    0
}

/// OCCT BSplCLib::FactorBandedMatrix (BSplCLib_2.cxx L385-...).
pub fn factor_banded_matrix(
    matrix: &mut BandMatrix,
    upper_band_width: i32,
    lower_band_width: i32,
    pivot_index_problem: &mut i32,
) -> i32 {
    let a_band_width = upper_band_width + lower_band_width + 1;
    *pivot_index_problem = 0;

    let l_row = 1i32;
    let n_cols = matrix.n_cols;
    for i in (l_row + 1)..=matrix.n_rows {
        let min_index = (lower_band_width - i + 2).max(1);

        for j in min_index..=lower_band_width {
            let index = i - lower_band_width + j - 1;
            // OCCT aRowIdx[LowerBandWidth] — 0-based column offset.
            let pivot = matrix.get(index, lower_band_width + 1);
            if pivot.abs() <= real_small() {
                *pivot_index_problem = index;
                return 1;
            }

            let inverse = -1.0 / pivot;
            // OCCT aRowI[j - 1] / aRowIdx[k + i - anIndex - 1] are 0-based
            // column offsets: BandMatrix column = offset + 1.
            let v = matrix.get(i, j) * inverse;
            matrix.set(i, j, v);
            let max_index = a_band_width + index - i;
            for k in (j + 1)..=max_index {
                let a = matrix.get(i, k);
                let m = matrix.get(i, j);
                let b = matrix.get(index, k + i - index);
                matrix.set(i, k, a + m * b);
            }
        }
    }

    0
}

/// OCCT BSplCLib::SolveBandedSystem (BSplCLib.cxx L3205-3270) — in-place
/// solve of the factored band system against `poles` (dimension-strided).
pub fn solve_banded_system(
    matrix: &BandMatrix,
    upper_band_width: i32,
    lower_band_width: i32,
    array_dimension: usize,
    poles: &mut [f64],
) -> i32 {
    if matrix.n_cols != upper_band_width + lower_band_width + 1 {
        return 1;
    }

    let l_row = 1i32;
    let dim = array_dimension as i32;
    for ii in (l_row + 1)..=matrix.n_rows {
        let min_index = (ii - lower_band_width).max(l_row);
        for jj in min_index..ii {
            let coeff = matrix.get(ii, jj - ii + lower_band_width + 1);
            for kk in 0..dim {
                poles[((ii - 1) * dim + kk) as usize] +=
                    poles[((jj - 1) * dim + kk) as usize] * coeff;
            }
        }
    }

    for ii in (l_row..=matrix.n_rows).rev() {
        let max_index = (ii + upper_band_width).min(matrix.n_rows);
        for jj in (ii + 1..=max_index).rev() {
            let coeff = matrix.get(ii, jj - ii + lower_band_width + 1);
            for kk in 0..dim {
                poles[((ii - 1) * dim + kk) as usize] -=
                    poles[((jj - 1) * dim + kk) as usize] * coeff;
            }
        }

        // Fixing a bug PRO18577 to avoid division by zero.
        let divisor = matrix.get(ii, lower_band_width + 1);
        const THE_TOLERANCE: f64 = 1.0e-16;
        if divisor.abs() <= THE_TOLERANCE {
            return 1;
        }
        let inverse = 1.0 / divisor;
        for kk in 0..dim {
            poles[((ii - 1) * dim + kk) as usize] *= inverse;
        }
    }

    0
}

/// OCCT BSplCLib::Interpolate (BSplCLib.cxx L3353-3395) — the non-homogeneous
/// overload.  On success `poles` holds the interpolating control values.
pub fn interpolate(
    degree: usize,
    flat_knots: &[f64],
    parameters: &[f64],
    contact_order_array: &[i32],
    array_dimension: usize,
    poles: &mut [f64],
) -> i32 {
    let mut upper_band_width = 0i32;
    let mut lower_band_width = 0i32;
    let mut inversion_problem = 0i32;
    let mut interpolation_matrix =
        BandMatrix::new(parameters.len() as i32, 2 * degree as i32 + 1);
    let error_code = build_bsp_matrix(
        parameters,
        contact_order_array,
        flat_knots,
        degree,
        &mut interpolation_matrix,
        &mut upper_band_width,
        &mut lower_band_width,
    );
    if error_code != 0 {
        return error_code;
    }
    let error_code = factor_banded_matrix(
        &mut interpolation_matrix,
        upper_band_width,
        lower_band_width,
        &mut inversion_problem,
    );
    if error_code != 0 {
        return error_code;
    }
    let error_code = solve_banded_system(
        &interpolation_matrix,
        upper_band_width,
        lower_band_width,
        array_dimension,
        poles,
    );
    if error_code != 0 {
        return error_code;
    }
    inversion_problem
}

/// OCCT BSplCLib::Eval — non-rational flat-knots overload (BSplCLib.cxx
/// L3640-3790).  `extrap_mode` is the 2-entry array whose FIRST entry is
/// passed by reference in OCCT (adjacent [1] is read inside Eval).
pub fn eval_flat(
    parameter: f64,
    periodic_flag: bool,
    derivative_request: i32,
    extrap_mode: &mut [i32; 2],
    degree: usize,
    flat_knots: &[f64],
    array_dimension: usize,
    poles: &[f64],
    results: &mut [f64],
) {
    let degree = degree as i32;
    let order = degree + 1;
    let mut local_request = derivative_request;
    let mut local_parameter = parameter;
    let mut extrapolating_flag = [0i32; 2];

    let k_upper = flat_knots.len() as i32;
    if periodic_flag {
        let period = at(flat_knots, k_upper - 1) - at(flat_knots, 2);
        while local_parameter > at(flat_knots, k_upper - 1) {
            local_parameter -= period;
        }
        while local_parameter < at(flat_knots, 2) {
            local_parameter += period;
        }
    }
    if parameter < at(flat_knots, 2)
        && local_request < extrap_mode[0]
        && extrap_mode[0] < degree
    {
        local_request = extrap_mode[0];
        local_parameter = at(flat_knots, 2);
        extrapolating_flag[0] = 1;
    }
    if parameter > at(flat_knots, k_upper - 1)
        && local_request < extrap_mode[1]
        && extrap_mode[1] < degree
    {
        local_request = extrap_mode[1];
        local_parameter = at(flat_knots, k_upper - 1);
        extrapolating_flag[1] = 1;
    }
    let delta = parameter - local_parameter;
    if local_request >= order {
        local_request = degree;
    }

    let modulus = if periodic_flag {
        (flat_knots.len() as i32) - degree - 1
    } else {
        (flat_knots.len() as i32) - degree
    };

    let mut first_non_zero = 0i32;
    let mut bspline_basis = BandMatrix::new(local_request + 1, order);
    let error_code = eval_bspline_basis(
        local_request,
        order,
        flat_knots,
        local_parameter,
        &mut first_non_zero,
        &mut bspline_basis,
        periodic_flag,
    );
    if error_code != 0 {
        return;
    }

    if extrapolating_flag[0] == 0 && extrapolating_flag[1] == 0 {
        let mut index = 0usize;
        for ii in 1..=(local_request + 1) {
            let mut index1 = first_non_zero;
            for kk in 0..array_dimension as i32 {
                results[index + kk as usize] = 0.0;
            }
            for jj in 1..=order {
                let b = bspline_basis.get(ii, jj);
                for kk in 0..array_dimension as i32 {
                    results[index + kk as usize] +=
                        poles[((index1 - 1) * array_dimension as i32 + kk) as usize] * b;
                }
                index1 = index1 % modulus;
                index1 += 1;
            }
            index += array_dimension;
        }
    } else {
        // Taylor expansion branch.
        let mut new_request = derivative_request;
        if new_request > degree {
            new_request = degree;
        }
        let n = ((local_request + 1) * array_dimension as i32) as usize;
        let mut local_real_array = vec![0.0f64; n];
        let mut index = 0usize;
        let mut inverse = 1.0f64;
        for ii in 1..=(local_request + 1) {
            let mut index1 = first_non_zero;
            for kk in 0..array_dimension as i32 {
                local_real_array[index + kk as usize] = 0.0;
            }
            for jj in 1..=order {
                let b = bspline_basis.get(ii, jj);
                for kk in 0..array_dimension as i32 {
                    local_real_array[index + kk as usize] +=
                        poles[((index1 - 1) * array_dimension as i32 + kk) as usize] * b;
                }
                index1 = index1 % modulus;
                index1 += 1;
            }
            for kk in 0..array_dimension as i32 {
                local_real_array[index + kk as usize] *= inverse;
            }
            index += array_dimension;
            inverse /= ii as f64;
        }
        eval_polynomial_flat(
            delta,
            new_request,
            degree,
            array_dimension as i32,
            &local_real_array,
            results,
        );
    }
}

/// OCCT BSplCLib::Eval — homogeneous flat-knots overload (BSplCLib.cxx
/// L3456-3638).  `weights` are the curve weights (unit weights for
/// non-rational curves, matching OCCT 8.0 semantics).
#[allow(clippy::too_many_arguments)]
pub fn eval_homogeneous(
    parameter: f64,
    periodic_flag: bool,
    derivative_request: i32,
    extrap_mode: &mut [i32; 2],
    degree: usize,
    flat_knots: &[f64],
    array_dimension: usize,
    poles: &[f64],
    weights: &[f64],
    poles_results: &mut [f64],
    weights_results: &mut [f64],
) {
    let degree = degree as i32;
    let order = degree + 1;
    let mut local_request = derivative_request;
    let mut local_parameter = parameter;
    let mut extrapolating_flag = [0i32; 2];

    let k_upper = flat_knots.len() as i32;
    if periodic_flag {
        let period = at(flat_knots, k_upper - 1) - at(flat_knots, 2);
        while local_parameter > at(flat_knots, k_upper - 1) {
            local_parameter -= period;
        }
        while local_parameter < at(flat_knots, 2) {
            local_parameter += period;
        }
    }
    if parameter < at(flat_knots, 2)
        && local_request < extrap_mode[0]
        && extrap_mode[0] < degree
    {
        local_request = extrap_mode[0];
        local_parameter = at(flat_knots, 2);
        extrapolating_flag[0] = 1;
    }
    if parameter > at(flat_knots, k_upper - 1)
        && local_request < extrap_mode[1]
        && extrap_mode[1] < degree
    {
        local_request = extrap_mode[1];
        local_parameter = at(flat_knots, k_upper - 1);
        extrapolating_flag[1] = 1;
    }
    let delta = parameter - local_parameter;
    if local_request >= order {
        local_request = degree;
    }

    let modulus = if periodic_flag {
        (flat_knots.len() as i32) - degree - 1
    } else {
        (flat_knots.len() as i32) - degree
    };

    let mut first_non_zero = 0i32;
    let mut bspline_basis = BandMatrix::new(local_request + 1, order);
    let error_code = eval_bspline_basis(
        local_request,
        order,
        flat_knots,
        local_parameter,
        &mut first_non_zero,
        &mut bspline_basis,
        periodic_flag,
    );
    if error_code != 0 {
        return;
    }

    if extrapolating_flag[0] == 0 && extrapolating_flag[1] == 0 {
        let mut index = 0usize;
        let mut index2 = 0usize;
        for ii in 1..=(local_request + 1) {
            let mut index1 = first_non_zero;
            for kk in 0..array_dimension as i32 {
                poles_results[index + kk as usize] = 0.0;
            }
            weights_results[index2] = 0.0;
            for jj in 1..=order {
                let b = bspline_basis.get(ii, jj);
                for kk in 0..array_dimension as i32 {
                    poles_results[index + kk as usize] += poles[((index1 - 1)
                        * array_dimension as i32
                        + kk) as usize]
                        * weights[(index1 - 1) as usize]
                        * b;
                }
                weights_results[index2] += weights[(index1 - 1) as usize] * b;
                index1 = index1 % modulus;
                index1 += 1;
            }
            index += array_dimension;
            index2 += 1;
        }
    } else {
        // Taylor expansion branch (poles, then weights).
        let mut new_request = derivative_request;
        if new_request > degree {
            new_request = degree;
        }
        let n = ((local_request + 1) * array_dimension as i32) as usize;
        let mut local_real_array = vec![0.0f64; n];
        let mut index = 0usize;
        let mut inverse = 1.0f64;
        for ii in 1..=(local_request + 1) {
            let mut index1 = first_non_zero;
            for kk in 0..array_dimension as i32 {
                local_real_array[index + kk as usize] = 0.0;
            }
            for jj in 1..=order {
                let b = bspline_basis.get(ii, jj);
                for kk in 0..array_dimension as i32 {
                    local_real_array[index + kk as usize] += poles[((index1 - 1)
                        * array_dimension as i32
                        + kk) as usize]
                        * weights[(index1 - 1) as usize]
                        * b;
                }
                index1 = index1 % modulus;
                index1 += 1;
            }
            for kk in 0..array_dimension as i32 {
                local_real_array[index + kk as usize] *= inverse;
            }
            index += array_dimension;
            inverse /= ii as f64;
        }
        eval_polynomial_flat(
            delta,
            new_request,
            degree,
            array_dimension as i32,
            &local_real_array,
            poles_results,
        );
        let mut index = 0usize;
        let mut inverse = 1.0f64;
        for ii in 1..=(local_request + 1) {
            let mut index1 = first_non_zero;
            local_real_array[index] = 0.0;
            for jj in 1..=order {
                let b = bspline_basis.get(ii, jj);
                local_real_array[index] += weights[(index1 - 1) as usize] * b;
                index1 = index1 % modulus;
                index1 += 1;
            }
            local_real_array[index] *= inverse;
            index += 1;
            inverse /= ii as f64;
        }
        eval_polynomial_flat(delta, new_request, degree, 1, &local_real_array, weights_results);
    }
}

/// OCCT PLib::RationalDerivatives (PLib.cxx) — converts homogeneous
/// derivatives into rational (divided) derivatives.  In the helix pipeline
/// (the only consumer) OCCT passes `RationalDerivates == PolesDerivates`
/// (one aliased buffer); this port therefore implements the in-place form:
/// `poles_and_rational` holds the homogeneous pole derivatives on entry and
/// the rational derivatives on exit.
pub fn rational_derivatives_inplace(
    derivative_request: i32,
    dimension: usize,
    poles_and_rational: &mut [f64],
    weights_derivates: &mut [f64],
) {
    let de_request1 = derivative_request + 1;
    let mut binomial_array = vec![1.0f64; de_request1 as usize];
    let inverse = 1.0 / weights_derivates[0];

    if dimension == 3 {
        let mut index = 0usize;
        let mut index2 = -6i32;
        for ii in 0..de_request1 {
            index2 += 3;
            let mut index1 = index2;
            // RationalArray[Index] = PolesArray[Index] — same buffer, no-op.
            index += 2;
            index -= 2;
            index += 2;
            index -= 2;
            for jj in (0..ii).rev() {
                let factor = binomial_array[jj as usize] * weights_derivates[(ii - jj) as usize];
                poles_and_rational[index] -= factor * poles_and_rational[index1 as usize];
                index += 1;
                index1 += 1;
                poles_and_rational[index] -= factor * poles_and_rational[index1 as usize];
                index += 1;
                index1 += 1;
                poles_and_rational[index] -= factor * poles_and_rational[index1 as usize];
                index -= 2;
                index1 -= 5;
            }
            for jj in (1..=ii).rev() {
                binomial_array[jj as usize] += binomial_array[(jj - 1) as usize];
            }
            poles_and_rational[index] *= inverse;
            index += 1;
            poles_and_rational[index] *= inverse;
            index += 1;
            poles_and_rational[index] *= inverse;
            index += 1;
        }
    } else {
        let dimension2 = (dimension << 1) as i32;
        let mut index = 0usize;
        let mut index2 = -dimension2;
        for ii in 0..de_request1 {
            index2 += dimension as i32;
            let mut index1 = index2;
            // RationalArray[Index] = PolesArray[Index] — same buffer.
            index += dimension;
            index -= dimension;
            for jj in (0..ii).rev() {
                let factor = binomial_array[jj as usize] * weights_derivates[(ii - jj) as usize];
                for _kk in 0..dimension {
                    poles_and_rational[index] -= factor * poles_and_rational[index1 as usize];
                    index += 1;
                    index1 += 1;
                }
                index -= dimension;
                index1 -= dimension2;
            }
            for jj in (1..=ii).rev() {
                binomial_array[jj as usize] += binomial_array[(jj - 1) as usize];
            }
            for _kk in 0..dimension {
                poles_and_rational[index] *= inverse;
                index += 1;
            }
        }
    }
}

/// OCCT BSplCLib::MovePointAndTangent (BSplCLib_2.cxx L567-864) — moves the
/// curve so that it passes through `delta` (offset from the current point)
/// with derivative `delta_derivatives` at U, disturbing only the poles in
/// [starting_condition, ending_condition] range.
#[allow(clippy::too_many_arguments)]
pub fn move_point_and_tangent(
    u: f64,
    array_dimension: usize,
    delta: &[f64],
    delta_derivatives: &[f64],
    tolerance: f64,
    degree: usize,
    starting_condition: i32,
    ending_condition: i32,
    poles: &[f64],
    weights: Option<&[f64]>,
    flat_knots: &[f64],
    new_poles: &mut [f64],
) -> i32 {
    let degree = degree as i32;
    let order = degree + 1;
    let num_knots = flat_knots.len() as i32;
    let num_poles = num_knots - order;
    let conditions = starting_condition + ending_condition + 4;

    if !(starting_condition >= -1 && starting_condition <= degree)
        || !(ending_condition >= -1 && ending_condition <= degree)
        || conditions > num_poles
    {
        return 2;
    }

    // check the parameter is within bounds.
    let mut start_index = 1 + degree;
    let mut end_index = num_knots - degree;
    let mut conditions_ok = true;
    if starting_condition == -1 {
        conditions_ok &= at(flat_knots, start_index) <= u;
    } else {
        conditions_ok &= at(flat_knots, start_index) + tolerance < u;
    }
    if ending_condition == -1 {
        conditions_ok &= at(flat_knots, end_index) >= u;
    } else {
        conditions_ok &= at(flat_knots, end_index) - tolerance > u;
    }
    if !conditions_ok {
        return 1;
    }

    // build 2 auxiliary functions.
    let mut schoenberg_points = vec![0.0f64; num_poles as usize];
    let mut first_function = vec![0.0f64; num_poles as usize];
    let mut second_function = vec![0.0f64; num_poles as usize];

    build_schoenberg_points(degree as usize, flat_knots, &mut schoenberg_points);
    let start_index = starting_condition + 2;
    let end_index = num_poles - ending_condition - 1;
    let mut index = 0i32;
    let mut new_parameter = 0.0f64;
    locate_parameter_main(
        &schoenberg_points,
        u,
        false,
        start_index,
        end_index,
        &mut index,
        &mut new_parameter,
        0.0,
        1.0,
    );

    let other_index;
    if index == start_index {
        other_index = index + 1;
    } else if index == end_index {
        other_index = index - 1;
    } else if u - at(&schoenberg_points, index) < at(&schoenberg_points, index + 1) - u {
        other_index = index - 1;
    } else {
        other_index = index + 1;
    }
    let type_: i32 = 3;

    let start_num_poles = starting_condition + 2;
    let end_num_poles = num_poles - ending_condition - 1;
    let start_value;
    if start_num_poles == 1 {
        let v = at(&schoenberg_points, num_poles) - at(&schoenberg_points, 1);
        start_value = at(&schoenberg_points, 1) - v;
    } else {
        start_value = at(&schoenberg_points, start_num_poles - 1);
    }
    let end_value;
    if end_num_poles == num_poles {
        let v = at(&schoenberg_points, num_poles) - at(&schoenberg_points, 1);
        end_value = at(&schoenberg_points, num_poles) + v;
    } else {
        end_value = at(&schoenberg_points, end_num_poles + 1);
    }

    for ii in 1..start_num_poles {
        first_function[(ii - 1) as usize] = 0.0;
        second_function[(ii - 1) as usize] = 0.0;
    }
    for ii in (end_num_poles + 1)..=num_poles {
        first_function[(ii - 1) as usize] = 0.0;
        second_function[(ii - 1) as usize] = 0.0;
    }
    let mut divide = 1.0 / (at(&schoenberg_points, index) - start_value);
    for ii in start_num_poles..=index {
        let mut value = at(&schoenberg_points, ii) - start_value;
        value *= divide;
        first_function[(ii - 1) as usize] = 1.0;
        for _jj in 0..type_ {
            first_function[(ii - 1) as usize] *= value;
        }
    }
    divide = 1.0 / (end_value - at(&schoenberg_points, index));
    for ii in index..=end_num_poles {
        let mut value = end_value - at(&schoenberg_points, ii);
        value *= divide;
        first_function[(ii - 1) as usize] = 1.0;
        for _jj in 0..type_ {
            first_function[(ii - 1) as usize] *= value;
        }
    }
    divide = 1.0 / (at(&schoenberg_points, other_index) - start_value);
    for ii in start_num_poles..=other_index {
        let mut value = at(&schoenberg_points, ii) - start_value;
        value *= divide;
        second_function[(ii - 1) as usize] = 1.0;
        for _jj in 0..type_ {
            second_function[(ii - 1) as usize] *= value;
        }
    }
    divide = 1.0 / (end_value - at(&schoenberg_points, other_index));
    for ii in other_index..=end_num_poles {
        let mut value = end_value - at(&schoenberg_points, ii);
        value *= divide;
        second_function[(ii - 1) as usize] = 1.0;
        for _jj in 0..type_ {
            second_function[(ii - 1) as usize] *= value;
        }
    }

    // compute the point and derivatives of both functions.
    let mut results = [[0.0f64; 2]; 2];
    let mut weights_results = [[0.0f64; 2]; 2];
    let derivative_request = 1;
    let dimension = 1usize;
    let mut extrap_mode = [degree, degree];

    match weights {
        Some(w) => {
            // evaluate in homogenised form.
            let mut r0 = [0.0f64; 2];
            let mut w0 = [0.0f64; 2];
            eval_homogeneous(
                u,
                false,
                derivative_request,
                &mut extrap_mode,
                degree as usize,
                flat_knots,
                dimension,
                &first_function,
                w,
                &mut r0,
                &mut w0,
            );
            let mut r1 = [0.0f64; 2];
            let mut w1 = [0.0f64; 2];
            eval_homogeneous(
                u,
                false,
                derivative_request,
                &mut extrap_mode,
                degree as usize,
                flat_knots,
                dimension,
                &second_function,
                w,
                &mut r1,
                &mut w1,
            );
            // compute the rational derivatives values.
            rational_derivatives_inplace(1, 1, &mut r0, &mut w0);
            rational_derivatives_inplace(1, 1, &mut r1, &mut w1);
            results[0] = r0;
            results[1] = r1;
        }
        None => {
            let mut r0 = [0.0f64; 2];
            eval_flat(
                u,
                false,
                1,
                &mut extrap_mode,
                degree as usize,
                flat_knots,
                1,
                &first_function,
                &mut r0,
            );
            let mut r1 = [0.0f64; 2];
            eval_flat(
                u,
                false,
                1,
                &mut extrap_mode,
                degree as usize,
                flat_knots,
                1,
                &second_function,
                &mut r1,
            );
            results[0] = r0;
            results[1] = r1;
        }
    }

    let m00 = results[0][0];
    let m01 = results[0][1];
    let m10 = results[1][0];
    let m11 = results[1][1];
    // a_matrix = inverse of [[results[0][0], results[0][1]], [results[1][0],
    // results[1][1]]]: inv = [[m11, -m01], [-m10, m00]] / det.
    let det = m00 * m11 - m01 * m10;
    if det == 0.0 {
        return 1;
    }
    let a00 = m11 / det;
    let a01 = -m01 / det;
    let a10 = -m10 / det;
    let a11 = m00 / det;

    let mut the_a_vector = vec![0.0f64; array_dimension];
    let mut the_b_vector = vec![0.0f64; array_dimension];
    for ii in 0..array_dimension {
        the_a_vector[ii] = a00 * delta[ii] + a10 * delta_derivatives[ii];
        the_b_vector[ii] = a01 * delta[ii] + a11 * delta_derivatives[ii];
    }

    let mut index = 0usize;
    for ii in 0..num_poles as usize {
        for jj in 0..array_dimension {
            new_poles[index] = poles[index];
            new_poles[index] += first_function[ii] * the_a_vector[jj];
            new_poles[index] += second_function[ii] * the_b_vector[jj];
            index += 1;
        }
    }

    0
}

// ---------------------------------------------------------------------------
// Knot insertion and degree elevation (BSplCLib.cxx L1555-2560)
// ---------------------------------------------------------------------------

/// OCCT BSplCLib::BuildKnots — builds the 2*Degree local knot array around
/// index `index`.  `mults` may be None (flat knots).
fn build_knots_local(
    degree: usize,
    index: i32,
    periodic: bool,
    knots: &[f64],
    mults: Option<&[i32]>,
    knot: &mut [f64],
) {
    let degree_i = degree as i32;
    match mults {
        None => {
            // The unrolled switch in OCCT copies knot[Index-Degree+i],
            // i = 0..2*Degree — a plain copy loop.
            let mut j = index - degree_i;
            for i in 0..(2 * degree) {
                j += 1;
                knot[i] = at(knots, j);
            }
        }
        Some(mults) => {
            let deg1 = degree_i - 1;
            let k_upper = knots.len() as i32;
            let m_upper = mults.len() as i32;
            let mut dknot = 0.0f64;
            let mut ilow = index;
            let mut mlow = 0;
            let mut iupp = index + 1;
            let mut mupp = 0;
            let mut loffset = 0.0f64;
            let mut uoffset = 0.0f64;
            let mut getlow = true;
            let mut getupp = true;
            if periodic {
                dknot = at(knots, k_upper) - at(knots, 1);
                if iupp > m_upper {
                    iupp = 2;
                    uoffset = dknot;
                }
            }
            // Find the knots around Index.
            for i in 0..degree {
                if getlow {
                    mlow += 1;
                    if mlow > ati(mults, ilow) {
                        mlow = 1;
                        ilow -= 1;
                        getlow = ilow >= 1;
                        if periodic && !getlow {
                            ilow = m_upper - 1;
                            loffset = dknot;
                            getlow = true;
                        }
                    }
                    if getlow {
                        knot[(deg1 - i as i32) as usize] = at(knots, ilow) - loffset;
                    }
                }
                if getupp {
                    mupp += 1;
                    if mupp > ati(mults, iupp) {
                        mupp = 1;
                        iupp += 1;
                        getupp = iupp <= m_upper;
                        if periodic && !getupp {
                            iupp = 2;
                            uoffset = dknot;
                            getupp = true;
                        }
                    }
                    if getupp {
                        knot[(degree_i + i as i32) as usize] = at(knots, iupp) + uoffset;
                    }
                }
            }
        }
    }
}

/// OCCT BSplCLib::BuildBoor.
fn build_boor(index: i32, length: i32, dimension: usize, poles: &[f64], lp: &mut [f64]) {
    let dimension = dimension as i32;
    let mut ip = index * dimension;
    let poles_len = poles.len() as i32;
    let mut ptr = 0usize;
    for _i in 0..=length {
        for k in 0..dimension {
            lp[ptr + k as usize] = poles[(ip - 1) as usize];
            ip += 1;
            if ip > poles_len {
                ip = 1;
            }
        }
        ptr += (2 * dimension) as usize;
    }
}

/// OCCT BSplCLib::BoorIndex.
fn boor_index(index: i32, length: i32, depth: i32) -> i32 {
    if index <= depth {
        return index;
    }
    if index <= length {
        return 2 * index - depth;
    }
    length + index - depth
}

/// OCCT BSplCLib::GetPole.
fn get_pole(
    index: i32,
    length: i32,
    depth: i32,
    dimension: usize,
    lp: &[f64],
    position: &mut i32,
    pole: &mut [f64],
) {
    let dimension = dimension as i32;
    let base = boor_index(index, length, depth) * dimension;
    for k in 0..dimension {
        set_at(pole, *position, lp[(base + k) as usize]);
        *position += 1;
        if *position > pole.len() as i32 {
            *position = 1;
        }
    }
}

/// OCCT BSplCLib::BoorScheme.
fn boor_scheme(
    u: f64,
    degree: usize,
    knots: &[f64],
    dimension: usize,
    poles: &mut [f64],
    depth: i32,
    length: i32,
) {
    let dimension = dimension as i32;
    // OCCT walks raw pointers: firstpole = &Poles - 2*Dimension, then
    // firstpole += Dimension per step and pole += 2*Dimension per i.  The
    // pointer offset o addresses the 0-based slot o (element Lower + o), so
    // the 0-based slot of `pole[k]` is pole_off + k.
    let dim = dimension;
    let mut firstpole: i32 = -2 * dim;
    for _step in 0..depth {
        firstpole += dim;
        let mut pole = firstpole;
        for i in _step..length {
            pole += 2 * dim;
            let x = (at(knots, i + degree as i32 - _step) - u)
                / (at(knots, i + degree as i32 - _step) - at(knots, i));
            let y = 1.0 - x;
            for k in 0..dim {
                let prev = poles[(pole + k - dim) as usize];
                let next = poles[(pole + k + dim) as usize];
                poles[(pole + k) as usize] = x * prev + y * next;
            }
        }
    }
}

/// OCCT static Copy (BSplCLib.cxx L2028-2054).
fn copy_poles(
    nb_poles: i32,
    old_first: &mut i32,
    old_poles: &[f64],
    new_first: &mut i32,
    new_poles: &mut [f64],
) {
    let old_lower = 1i32;
    let old_upper = old_poles.len() as i32;
    let new_lower = 1i32;
    let new_upper = new_poles.len() as i32;

    *old_first = old_lower + (*old_first - old_lower).rem_euclid(old_upper - old_lower + 1);
    *new_first = new_lower + (*new_first - new_lower).rem_euclid(new_upper - new_lower + 1);

    for _i in 1..=nb_poles {
        set_at(new_poles, *new_first, at(old_poles, *old_first));
        *old_first += 1;
        if *old_first > old_upper {
            *old_first = old_lower;
        }
        *new_first += 1;
        if *new_first > new_upper {
            *new_first = new_lower;
        }
    }
}

/// OCCT BSplCLib::InsertKnots (BSplCLib.cxx L2063-2355).  `add_mults == None`
/// means `addflat` (insert each added knot once).
#[allow(clippy::too_many_arguments)]
pub fn insert_knots(
    degree: usize,
    periodic: bool,
    dimension: usize,
    poles: &[f64],
    knots: &[f64],
    mults: &[i32],
    add_knots: &[f64],
    add_mults: Option<&[i32]>,
    new_poles: &mut [f64],
    new_knots: &mut [f64],
    new_mults: &mut [i32],
    tolerance: f64,
    add: bool,
) {
    let addflat = add_mults.is_none();
    let degree_i = degree as i32;
    let dim = dimension as i32;

    let mut knots_local = vec![0.0f64; 2 * degree];
    let mut poles_local = vec![0.0f64; ((2 * degree + 1) * dimension)];

    let knots_upper = knots.len() as i32;
    let mut curk = 0i32; // Knots.Lower() - 1
    let mut curnk = 0i32; // NewKnots.Lower() - 1
    let mut curp = 1i32;
    let mut curnp = 1i32;

    let mut index;
    if periodic {
        index = -ati(mults, 1);
    } else {
        index = -degree_i - 1;
    }

    let mut firstmult = 0i32;

    for kn in 1..=add_knots.len() as i32 {
        let u = at(add_knots, kn);
        let mut eps = tolerance.max(epsilon() * u.abs());

        // find the position in the old knots and copy to the new knots.
        while curk < knots_upper && at(knots, curk + 1) - u <= eps {
            curk += 1;
            curnk += 1;
            set_at(new_knots, curnk, at(knots, curk));
            let m = ati(mults, curk);
            set_ati(new_mults, curnk, m);
            index += m;
        }

        // Slice the knots and the mults to the current size of the new curve.
        let slice_len = (curnk + knots_upper - curk) as usize;
        // nknots/nmults are views into NewKnots/NewMults with length
        // slice_len (1-based indexing preserved via at/set_at).
        // Copy enough knots to compute the insertion schema.
        let mut k = curk;
        let mut i = curnk;
        let mut mult = 0i32;
        while (mult as usize) < degree && k < knots_upper {
            k += 1;
            i += 1;
            set_at(new_knots, i, at(knots, k));
            let m = ati(mults, k);
            set_ati(new_mults, i, m);
            mult += m;
        }

        // copy knots at the end for periodic curve.
        if periodic {
            let mut mult = 0i32;
            let mut k = knots_upper;
            let mut i = slice_len as i32;
            while (mult as usize) < degree && i > curnk {
                set_at(new_knots, i, at(knots, k));
                let m = ati(mults, k);
                set_ati(new_mults, i, m);
                mult += m;
                k -= 1;
                i -= 1;
            }
            set_ati(new_mults, slice_len as i32, ati(new_mults, 1));
        }

        // do the boor scheme on the new curve to insert the new knot.
        let sameknot = (u - at(new_knots, curnk.max(1))).abs() <= eps;

        let length;
        if sameknot {
            length = (degree_i - ati(new_mults, curnk.max(1))).max(0);
        } else {
            length = degree_i;
        }

        let mut depth;
        if addflat {
            depth = 1;
        } else {
            depth = degree_i.min(ati(add_mults.unwrap(), kn));
        }

        if sameknot {
            if add {
                if ati(new_mults, curnk.max(1)) + depth > degree_i {
                    depth = degree_i - ati(new_mults, curnk.max(1));
                }
            } else {
                depth = (depth - ati(new_mults, curnk.max(1))).max(0);
            }

            if periodic {
                // on periodic curve the first and last knot are delayed to the end.
                if curk == 1 || curk == knots_upper {
                    if firstmult == 0 {
                        firstmult += depth;
                    }
                    depth = 0;
                }
            }
        }
        if depth <= 0 {
            continue;
        }

        build_knots_local(
            degree,
            curnk.max(1),
            periodic,
            new_knots,
            Some(new_mults),
            &mut knots_local,
        );

        // copy the poles.
        let mut need = 1 + (index + length + 1) * dim - curnp;
        need = need.min(poles.len() as i32 - curp + 1);

        let mut p = curp;
        let mut np = curnp;
        copy_poles(need, &mut p, poles, &mut np, new_poles);
        curp = p;
        curnp = np;

        // slice the poles to the current number of poles in case of periodic.
        // npoles view: NewPoles[1..=curnp-1].
        build_boor(
            index,
            length,
            dimension,
            &new_poles[..(curnp - 1).max(0) as usize],
            &mut poles_local,
        );
        boor_scheme(u, degree, &knots_local, dimension, &mut poles_local, depth, length);

        // copy the new poles.
        curnp += depth * dim;
        let mut np = 1 + (index + 1) * dim;
        for i in 1..=(length + depth) {
            get_pole(
                i,
                length,
                depth,
                dimension,
                &poles_local,
                &mut np,
                new_poles,
            );
        }

        // insert the knot.
        index += depth;
        if sameknot {
            let m = ati(new_mults, curnk.max(1));
            set_ati(new_mults, curnk.max(1), m + depth);
        } else {
            curnk += 1;
            set_at(new_knots, curnk, u);
            set_ati(new_mults, curnk, depth);
        }
        let _ = eps;
    }

    // copy the last poles and knots.
    let mut p = curp;
    let mut np = curnp;
    let need = poles.len() as i32 - curp + 1;
    copy_poles(need, &mut p, poles, &mut np, new_poles);
    curnp = np;

    while curk < knots_upper {
        curk += 1;
        curnk += 1;
        set_at(new_knots, curnk, at(knots, curk));
        set_ati(new_mults, curnk, ati(mults, curk));
    }

    // process the first-last knot on periodic curves.
    if firstmult > 0 {
        curnk = 1;
        if ati(new_mults, 1) + firstmult > degree_i {
            firstmult = degree_i - ati(new_mults, 1);
        }
        if firstmult > 0 {
            let length = degree_i - ati(new_mults, 1);
            let depth = firstmult;

            build_knots_local(degree, 1, periodic, new_knots, Some(new_mults), &mut knots_local);
            let npoles_upper = new_poles.len() as i32 - depth * dim;
            build_boor(
                0,
                length,
                dimension,
                &new_poles[..npoles_upper.max(0) as usize],
                &mut poles_local,
            );
            boor_scheme(
                at(new_knots, 1),
                degree,
                &knots_local,
                dimension,
                &mut poles_local,
                depth,
                length,
            );

            // copy the new poles but rotate them with depth.
            let mut np = 1i32;
            for i in depth..(length + depth) {
                get_pole(i, length, depth, dimension, &poles_local, &mut np, new_poles);
            }
            let mut np = new_poles.len() as i32 - depth * dim + 1;
            for i in 0..depth {
                get_pole(i, length, depth, dimension, &poles_local, &mut np, new_poles);
            }

            let m = ati(new_mults, 1);
            set_ati(new_mults, 1, m + depth);
            let m_last = ati(new_mults, new_mults.len() as i32);
            set_ati(new_mults, new_mults.len() as i32, m_last + depth);
        }
    }
}

/// OCCT BSplCLib::IncreaseDegreeCountKnots (BSplCLib.cxx).
pub fn increase_degree_count_knots(
    degree: usize,
    new_degree: usize,
    periodic: bool,
    mults: &[i32],
) -> usize {
    if periodic {
        return mults.len();
    }
    let degree_i = degree as i32;
    let new_degree_i = new_degree as i32;
    let f = first_uknot_index_mults(degree, mults) as i32;
    let l = last_uknot_index_mults(degree, mults) as i32;
    let step = (new_degree - degree) as i32;

    let mut removed = 0i32;
    let mut i = 1i32; // Mults.Lower()
    let mut m = degree_i + (f - i + 1) * step + 1;
    while m > new_degree_i + 1 {
        removed += 1;
        m -= ati(mults, i) + step;
        i += 1;
    }
    if m < new_degree_i + 1 {
        removed -= 1;
    }

    let mut i = mults.len() as i32; // Mults.Upper()
    let mut m = degree_i + (i - l + 1) * step + 1;
    while m > new_degree_i + 1 {
        removed += 1;
        m -= ati(mults, i) + step;
        i -= 1;
    }
    if m < new_degree_i + 1 {
        removed -= 1;
    }

    (mults.len() as i32 - removed) as usize
}


/// OCCT BSplCLib::IncreaseDegree (BSplCLib.cxx L2592-2560) — degree elevation
/// of a BSpline curve (Prautzsch averaging).
#[allow(clippy::too_many_arguments)]
pub fn increase_degree(
    degree: usize,
    new_degree: usize,
    periodic: bool,
    dimension: usize,
    poles: &[f64],
    knots: &[f64],
    mults: &[i32],
    new_poles: &mut [f64],
    new_knots: &mut [f64],
    new_mults: &mut [i32],
) {
    let dim = dimension as i32;
    let mut pf: i32 = 0;
    let mut pl: i32 = 0;

    let mut nbwknots = knots.len();
    let f = first_uknot_index_mults(degree, mults) as i32;
    let l = last_uknot_index_mults(degree, mults) as i32;

    if periodic {
        nbwknots += (f - 1) as usize;
        pf = -(degree as i32) - 1;
        for i in 1..=f {
            pf += ati(mults, i);
        }
        nbwknots += (mults.len() as i32 - l) as usize;
        pl = -(degree as i32) - 1;
        for i in l..=mults.len() as i32 {
            pl += ati(mults, i);
        }
    }

    // copy the knots and multiplicities.
    let mut wknots = vec![0.0f64; nbwknots];
    let mut wmults = vec![0i32; nbwknots];
    if !periodic {
        wknots.copy_from_slice(knots);
        wmults.copy_from_slice(mults);
    } else {
        let period = at(knots, knots.len() as i32) - at(knots, 1);
        let mut i = 0usize;
        let knots_upper = knots.len() as i32;
        // OCCT: for (k = l; k < Knots.Upper(); k++) — forward.
        for k in l..knots_upper {
            i += 1;
            wknots[i - 1] = at(knots, k) - period;
            wmults[i - 1] = ati(mults, k);
        }
        for k in 1..=knots_upper {
            i += 1;
            wknots[i - 1] = at(knots, k);
            wmults[i - 1] = ati(mults, k);
        }
        // OCCT: for (k = Knots.Lower() + 1; k <= f; k++).
        for k in 2..=f {
            i += 1;
            wknots[i - 1] = at(knots, k) + period;
            wmults[i - 1] = ati(mults, k);
        }
    }

    // set the first and last mults to Degree+1 and add null poles.
    pf += degree as i32 + 1 - wmults[0];
    wmults[0] = degree as i32 + 1;
    pl += degree as i32 + 1 - wmults[nbwknots - 1];
    wmults[nbwknots - 1] = degree as i32 + 1;

    // poles of the working curve.
    let mut nbwpoles = 0i32;
    for i in 0..nbwknots {
        nbwpoles += wmults[i];
    }
    nbwpoles -= degree as i32 + 1;

    let wp_max = (nbwpoles + (nbwknots as i32 - 1) * (new_degree - degree) as i32) as usize;
    let mut wpoles = vec![0.0f64; wp_max * dimension];

    for i in 0..(pf * dim) as usize {
        wpoles[i] = 0.0;
    }

    let mut k = 1usize; // Poles.Lower()
    for i in (pf * dim) as usize..((nbwpoles - pl) * dim) as usize {
        wpoles[i] = poles[k - 1];
        k += 1;
        if k > poles.len() {
            k = 1;
        }
    }
    for i in ((nbwpoles - pl) * dim) as usize..(nbwpoles * dim) as usize {
        wpoles[i] = 0.0;
    }

    // temporary arrays.
    let mut tempc1 = vec![0.0f64; wp_max * dimension];
    let mut tempc2 = vec![0.0f64; wp_max * dimension];
    let mut iknots = vec![0.0f64; nbwknots];
    let mut nwpoles = vec![0.0f64; wp_max * dimension];

    // loop on degree incrementation.
    let mut nbp;
    let mut nbwp = nbwpoles as usize;

    for cur_deg in degree..new_degree {
        nbp = nbwp;
        nbwp = nbp + nbwknots - 1;

        for idx in 0..nbwp * dimension {
            nwpoles[idx] = 0.0;
        }

        for step in 0..=cur_deg {
            if step != 0 {
                for i in 0..nbwknots {
                    wmults[i] -= 1;
                }
            }

            let mut offset: usize = 0;
            for i in 0..nbp {
                offset += 1;
                for k in 0..dimension {
                    tempc1[(offset - 1) * dimension + k] = wpoles[i * dimension + k];
                }
                if i % (cur_deg + 1) == step {
                    offset += 1;
                    for k in 0..dimension {
                        tempc1[(offset - 1) * dimension + k] = wpoles[i * dimension + k];
                    }
                }
            }

            // knot multiplicities increased / knots inserted.
            let mut stepmult = step as i32 + 1;
            let mut nbknots = 0i32;
            let mut smult = 0i32;
            for k in 0..nbwknots {
                smult += wmults[k];
                if smult >= stepmult {
                    stepmult += cur_deg as i32 + 1;
                    wmults[k] += 1;
                } else {
                    nbknots += 1;
                    iknots[(nbknots - 1) as usize] = wknots[k];
                }
            }

            if nbknots > 0 {
                let aknots = iknots[..nbknots as usize].to_vec();
                let mut ncurve = vec![0.0f64; nbwp * dimension];
                let mut n_knots = vec![0.0f64; nbwknots];
                // OCCT passes wmults as BOTH input mults and output NewMults
                // (in-place update); clone the input to satisfy borrows.
                let wmults_in = wmults.clone();
                insert_knots(
                    cur_deg + 1,
                    false,
                    dimension,
                    &tempc1[..offset * dimension],
                    &wknots,
                    &wmults_in,
                    &aknots,
                    None,
                    &mut ncurve,
                    &mut n_knots,
                    &mut wmults,
                    0.0,
                    false,
                );
                for i in 0..nbwp * dimension {
                    nwpoles[i] += ncurve[i];
                }
            } else {
                for i in 0..nbwp * dimension {
                    nwpoles[i] += tempc1[i];
                }
            }
        }

        // The result is the average.
        for i in 0..nbwp * dimension {
            wpoles[i] = nwpoles[i] / (cur_deg + 1) as f64;
        }
    }

    // copy the results.
    let firstknot = if periodic {
        (mults.len() as i32 - l) as usize
    } else {
        f as usize
    };

    let mut m = 0i32;
    for k in 0..firstknot {
        m += wmults[k];
    }

    let mut k = 1usize;
    pf = 0;
    while m > new_degree as i32 + 1 {
        k += 1;
        m -= wmults[k - 1];
        pf += wmults[k - 1];
    }
    if m < new_degree as i32 + 1 {
        k -= 1;
        wmults[k - 1] += m - new_degree as i32 - 1;
        pf += m - new_degree as i32 - 1;
    }
    if periodic {
        k = firstknot;
    }

    // copy knots.
    for i in 0..new_knots.len() {
        new_knots[i] = wknots[k - 1 + i];
        new_mults[i] = wmults[k - 1 + i];
    }

    // copy poles.
    let mut pf_idx = (pf as usize) * dimension;
    for i in 0..new_poles.len() {
        pf_idx += 1;
        new_poles[i] = wpoles[pf_idx - 1];
    }
    let _ = dim;
}

// ---------------------------------------------------------------------------
// Kernels for Segment / SetNotPeriodic / knot distribution analysis
// (BSplCLib.cxx), consumed by Geom2d_BSplineCurve operations.
// ---------------------------------------------------------------------------

/// OCCT BSplCLib::LocateParameter(Degree, Knots, Mults, U, IsPeriodic, FromK1,
/// ToK2, KnotIndex, NewU) (BSplCLib.cxx L168-185).  Degree and Mults are
/// unused by OCCT itself (kept for signature parity).
pub fn locate_parameter_knots_mults(
    _degree: usize,
    knots: &[f64],
    _mults: &[i32],
    u: f64,
    is_periodic: bool,
    from_k1: i32,
    to_k2: i32,
    knot_index: &mut i32,
    new_u: &mut f64,
) {
    let (uf, ul) = if is_periodic {
        (at(knots, 1), at(knots, knots.len() as i32))
    } else {
        (0.0f64, 1.0f64)
    };
    locate_parameter_main(knots, u, is_periodic, from_k1, to_k2, knot_index, new_u, uf, ul);
}

/// OCCT BSplCLib::PoleIndex (BSplCLib.cxx L1758-1779).
pub fn pole_index(degree: usize, index: i32, periodic: bool, mults: &[i32]) -> i32 {
    let mut pindex = 0i32;
    for i in 1..=index {
        pindex += ati(mults, i);
    }
    if periodic {
        pindex -= ati(mults, 1);
    } else {
        pindex -= degree as i32 + 1;
    }
    pindex
}

/// OCCT BSplCLib::PrepareInsertKnots (BSplCLib.cxx L1849-2024).  `add_mults
/// == None` means `addflat` (insert each added knot once).  On success
/// `nb_poles` / `nb_knots` hold the new curve sizes.
#[allow(clippy::too_many_arguments)]
pub fn prepare_insert_knots(
    degree: usize,
    periodic: bool,
    knots: &[f64],
    mults: &[i32],
    add_knots: &[f64],
    add_mults: Option<&[i32]>,
    nb_poles: &mut i32,
    nb_knots: &mut i32,
    tolerance: f64,
    add: bool,
) -> bool {
    let addflat = add_mults.is_none();
    let degree_i = degree as i32;
    let knots_upper = knots.len() as i32;
    let add_upper = add_knots.len() as i32;

    let (first, last) = if periodic {
        (1i32, knots_upper)
    } else {
        (
            first_uknot_index_mults(degree, mults),
            last_uknot_index_mults(degree, mults),
        )
    };
    let adelta_k1 = at(knots, first) - at(add_knots, 1);
    let adelta_k2 = at(add_knots, add_upper) - at(knots, last);
    if adelta_k1 > tolerance {
        return false;
    }
    if adelta_k2 > tolerance {
        return false;
    }

    let mut sigma = 0i32;
    let mut amult;
    *nb_knots = 0;
    let mut k = 0i32; // Knots.Lower() - 1
    let mut ak = 1i32; // AddKnots.Lower()

    if periodic && add_upper > 1 {
        // gka for case when segments was produced on full period only one knot
        // was added in the end of curve
        if adelta_k1.abs() <= GP_RESOLUTION && adelta_k2.abs() <= GP_RESOLUTION {
            ak += 1;
        }
    }

    let mut a_last_knot_mult = ati(mults, knots_upper);
    let mut au;
    let mut oldau = at(add_knots, ak);
    let mut eps;

    while ak <= add_upper {
        au = at(add_knots, ak);
        if au < oldau {
            return false;
        }
        oldau = au;

        eps = tolerance.max(epsilon() * au.abs());

        while k < knots_upper && at(knots, k + 1) - au <= eps {
            k += 1;
            *nb_knots += 1;
            sigma += ati(mults, k);
        }

        if addflat {
            amult = 1;
        } else {
            amult = ati(add_mults.unwrap(), ak).max(0);
        }

        while ak < add_upper && (au - at(add_knots, ak + 1)).abs() <= eps {
            ak += 1;
            if add {
                if addflat {
                    amult += 1;
                } else {
                    amult += ati(add_mults.unwrap(), ak).max(0);
                }
            }
        }

        if (au - at(knots, k)).abs() <= eps {
            // identic to existing knot
            let mult = ati(mults, k);
            if add {
                if mult + amult > degree_i {
                    amult = (degree_i - mult).max(0);
                }
                sigma += amult;
            } else if amult > mult {
                if amult > degree_i {
                    amult = degree_i;
                }
                if k == knots_upper && periodic {
                    a_last_knot_mult = amult.max(mult);
                    sigma += 2 * (a_last_knot_mult - mult);
                } else {
                    sigma += amult - mult;
                }
            }
        } else {
            // not identic to existing knot
            if amult > 0 {
                if amult > degree_i {
                    amult = degree_i;
                }
                *nb_knots += 1;
                sigma += amult;
            }
        }

        ak += 1;
    }

    // count the last knots
    while k < knots_upper {
        k += 1;
        *nb_knots += 1;
        sigma += ati(mults, k);
    }

    if periodic {
        // for periodic B-Spline the requirement is that multiplicities of the
        // first and last knots must be equal.
        *nb_poles = sigma - a_last_knot_mult;
    } else {
        *nb_poles = sigma - degree_i - 1;
    }

    true
}

/// OCCT BSplCLib::PrepareUnperiodize (BSplCLib.cxx L2967-3020).
pub fn prepare_unperiodize(degree: usize, mults: &[i32], nb_knots: &mut i32, nb_poles: &mut i32) {
    let degree_i = degree as i32;
    let l = mults.len() as i32;
    // initialize NbKnots and NbPoles
    *nb_knots = l;
    *nb_poles = -degree_i - 1;
    for i in 1..=l {
        *nb_poles += ati(mults, i);
    }

    let mut sigma;
    let mut k;
    // Add knots at the beginning of the curve to raise multiplicities
    // to Degree + 1.
    sigma = ati(mults, 1);
    k = l - 1;
    while sigma < degree_i + 1 {
        sigma += ati(mults, k);
        *nb_poles += ati(mults, k);
        k -= 1;
        *nb_knots += 1;
    }
    // We must add exactly until Degree + 1 -> suppress the excedent.
    if sigma > degree_i + 1 {
        *nb_poles -= sigma - degree_i - 1;
    }

    // Add knots at the end of the curve to raise multiplicities
    // to Degree + 1.
    sigma = ati(mults, l);
    k = 2;
    while sigma < degree_i + 1 {
        sigma += ati(mults, k);
        *nb_poles += ati(mults, k);
        k += 1;
        *nb_knots += 1;
    }
    if sigma > degree_i + 1 {
        *nb_poles -= sigma - degree_i - 1;
    }
}

/// OCCT BSplCLib::Unperiodize (BSplCLib.cxx L3024-3080).  The Dimension
/// parameter is ignored by OCCT (plain pole copy); `poles` / `new_poles` are
/// flat dimension-strided arrays sized by [`prepare_unperiodize`].
pub fn unperiodize(
    degree: usize,
    mults: &[i32],
    knots: &[f64],
    poles: &[f64],
    new_mults: &mut [i32],
    new_knots: &mut [f64],
    new_poles: &mut [f64],
) {
    let degree_i = degree as i32;
    let l = mults.len() as i32;
    let poles_len = poles.len() as i32;
    let mut index = 0i32;
    // evaluation of index : number of knots to insert before knot(1) to
    // raise sum of multiplicities to <Degree + 1>
    let mut sigma = ati(mults, 1);
    let mut k = l - 1;
    while sigma < degree_i + 1 {
        sigma += ati(mults, k);
        k -= 1;
        index += 1;
    }

    let period = at(knots, l) - at(knots, 1);

    // set the 'interior' knots;
    for k in 1..=l {
        set_at(new_knots, k + index, at(knots, k));
        set_ati(new_mults, k + index, ati(mults, k));
    }

    // set the 'starting' knots;
    for k in 1..=index {
        set_at(new_knots, k, at(new_knots, k + l - 1) - period);
        set_ati(new_mults, k, ati(new_mults, k + l - 1));
    }
    set_ati(new_mults, 1, ati(new_mults, 1) - (sigma - degree_i - 1));

    // set the 'ending' knots;
    sigma = ati(new_mults, index + l);
    let new_upper = new_knots.len() as i32;
    for k in (l + index + 1)..=new_upper {
        set_at(new_knots, k, at(new_knots, k - l + 1) + period);
        set_ati(new_mults, k, ati(new_mults, k - l + 1));
        sigma += ati(new_mults, k - l + 1);
    }
    let nm_upper = new_mults.len() as i32;
    set_ati(new_mults, nm_upper, ati(new_mults, nm_upper) - (sigma - degree_i - 1));

    for k in 1..=(new_poles.len() as i32) {
        set_at(new_poles, k, at(poles, (k - 1) % poles_len + 1));
    }
}

/// OCCT GeomAbs_BSplKnotDistribution.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GeomAbsKnotDistribution {
    NonUniform,
    Uniform,
    QuasiUniform,
    PiecewiseBezier,
    PiecewiseBezierAndPeriodic,
}

/// OCCT BSplCLib_KnotDistribution (internal KnotForm result).
#[derive(Clone, Copy, PartialEq, Eq)]
enum KnotFormDist {
    Uniform,
    NonUniform,
}

/// OCCT BSplCLib_MultDistribution (internal MultForm result).
#[derive(Clone, Copy, PartialEq, Eq)]
enum MultFormDist {
    Constant,
    NonConstant,
    QuasiConstant,
}

/// OCCT BSplCLib::KnotForm (BSplCLib.cxx L602-632).
fn knot_form(knots: &[f64], from_k1: i32, to_k2: i32) -> KnotFormDist {
    if from_k1 + 1 > knots.len() as i32 {
        return KnotFormDist::Uniform;
    }

    let mut a_ui = at(knots, from_k1).abs();
    let mut a_uj = at(knots, from_k1 + 1).abs();
    let mut a_du0 = (a_uj - a_ui).abs();
    let mut an_eps = epsilon() * a_ui + epsilon() * a_uj + epsilon() * a_du0;

    for i in (from_k1 + 1)..to_k2 {
        a_ui = at(knots, i).abs();
        a_uj = at(knots, i + 1).abs();
        let a_du1 = (a_uj - a_ui).abs();

        if (a_du1 - a_du0).abs() > an_eps {
            return KnotFormDist::NonUniform;
        }

        a_du0 = a_du1;
        an_eps = epsilon() * a_ui + epsilon() * a_uj + epsilon() * a_du0;
    }

    KnotFormDist::Uniform
}

/// OCCT BSplCLib::MultForm (BSplCLib.cxx L636-685).
fn mult_form(mults: &[i32], from_k1: i32, to_k2: i32) -> MultFormDist {
    let a_first = from_k1.min(to_k2);
    let a_last = from_k1.max(to_k2);

    if a_first + 1 > mults.len() as i32 {
        return MultFormDist::Constant;
    }

    let a_first_mult = ati(mults, a_first);
    let mut a_form = MultFormDist::Constant;
    let mut a_mult = ati(mults, a_first + 1);

    let mut i = a_first + 1;
    while i <= a_last && a_form != MultFormDist::NonConstant {
        if i == a_first + 1 {
            if a_mult != a_first_mult {
                a_form = MultFormDist::QuasiConstant;
            }
        } else if i == a_last {
            if matches!(a_form, MultFormDist::QuasiConstant) {
                if a_first_mult != ati(mults, i) {
                    a_form = MultFormDist::NonConstant;
                }
            } else if a_mult != ati(mults, i) {
                a_form = MultFormDist::NonConstant;
            }
        } else {
            if a_mult != ati(mults, i) {
                a_form = MultFormDist::NonConstant;
            }
            a_mult = ati(mults, i);
        }
        i += 1;
    }

    a_form
}

/// OCCT BSplCLib::KnotAnalysis (BSplCLib.cxx L692-755).
pub fn knot_analysis(
    degree: usize,
    periodic: bool,
    knots: &[f64],
    mults: &[i32],
    knot_form_out: &mut GeomAbsKnotDistribution,
    max_knot_mult: &mut i32,
) {
    let degree_i = degree as i32;
    *knot_form_out = GeomAbsKnotDistribution::NonUniform;

    let k_set = knot_form(knots, 1, knots.len() as i32);

    if matches!(k_set, KnotFormDist::Uniform) {
        match mult_form(mults, 1, mults.len() as i32) {
            MultFormDist::NonConstant => {}
            MultFormDist::Constant => {
                if knots.len() == 2 {
                    *knot_form_out = GeomAbsKnotDistribution::PiecewiseBezier;
                } else if ati(mults, 1) == 1 {
                    *knot_form_out = GeomAbsKnotDistribution::Uniform;
                }
            }
            MultFormDist::QuasiConstant => {
                if ati(mults, 1) == degree_i + 1 {
                    let m = ati(mults, 2);
                    if m == degree_i {
                        *knot_form_out = GeomAbsKnotDistribution::PiecewiseBezier;
                    } else if m == 1 {
                        *knot_form_out = GeomAbsKnotDistribution::QuasiUniform;
                    }
                }
            }
        }
    }

    let first_km = if periodic {
        1i32
    } else {
        first_uknot_index_mults(degree, mults)
    };
    let last_km = if periodic {
        knots.len() as i32
    } else {
        last_uknot_index_mults(degree, mults)
    };
    *max_knot_mult = 0;
    if last_km - first_km != 1 {
        for i in (first_km + 1)..last_km {
            let multi = ati(mults, i);
            *max_knot_mult = (*max_knot_mult).max(multi);
        }
    }
}
