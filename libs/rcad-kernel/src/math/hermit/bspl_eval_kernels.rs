//! Shared BSplCLib point/scalar evaluation kernels (OCCT TKMath/BSplCLib),
//! hosted inside the `hermit` submodule because this batch's edit surface on
//! `math/bspl_lib.rs` is fixed (FunctionReparameterise + periodic-arm helpers
//! only).  Both `Hermit` (scalar `BSplCLib::D1`) and
//! `Geom2d_BSplineCurve::EvalD0/EvalD1` (point `BSplCLib::D0/D1`) consume
//! these statics, exactly as their OCCT counterparts do.
//!
//! OCCT anchors:
//! - `BSplCLib::LocateParameter` (mults, default-range) BSplCLib.cxx L321-364
//! - `BSplCLib::BuildKnots`                       BSplCLib.cxx L1555-1740
//! - `BSplCLib::IsRational`                       BSplCLib.cxx L842-861
//! - `BSplCLib::BuildEval` (scalar)               BSplCLib_2.cxx L49-89
//! - `BSplCLib::PrepareEval` (scalar, file static) BSplCLib_2.cxx L93-138
//! - `BSplCLib::D0/D1` (scalar)                   BSplCLib_2.cxx L142-197
//! - `validateBSplineDegree`                      BSplCLib_CurveComputation.pxx L254-258
//! - `BSplCLib_BuildEval<Point,...,2>`            BSplCLib_CurveComputation.pxx L720-762
//! - `PrepareEval_T<Point,...,2>`                 BSplCLib_CurveComputation.pxx L776-830
//! - `BSplCLib_D0<Point,...,2>`                   BSplCLib_CurveComputation.pxx L843-882
//! - `BSplCLib_D1<Point,...,2>`                   BSplCLib_CurveComputation.pxx L896-935
//! - `BSplCLib::Bohm` (Dimension 1 / 2 / 3 cases) BSplCLib.cxx L1197-1401
//! - `PLib::RationalDerivative`                   PLib.cxx L274-558
//!
//! `BSplCLib::Eval` (in-place corner cutting, BSplCLib.cxx L865-1006) already
//! lives in [`crate::math::bspl_lib::bspl_clib_eval_inplace`] and is reused.
//!
//! Architecture note: OCCT `BSplCLib_DataContainer_T<N>` (pxx L263-268) is a
//! stack struct `poles[(MaxDegree+1)*(N+1)], knots[2*MaxDegree], ders[N*4]`
//! with `THE_MAX_DEGREE == BSplCLib::MaxDegree() == 25`; the Rust mirror below
//! keeps the identical layout with fixed-size arrays.  Pointer walks over the
//! local buffers use `isize` offsets exactly as OCCT moves raw pointers
//! (the walks may step one slot before the buffer between iterations without
//! dereferencing it); every dereference stays inside the buffer.

use crate::math::bspl_lib::{
    at, ati, bspl_clib_eval_inplace, first_uknot_index_mults, last_uknot_index_mults,
    locate_parameter_main, pole_index, BSPLIB_MAX_DEGREE,
};

/// OCCT `validateBSplineDegree` (pxx L254-258).
#[inline]
fn validate_bspline_degree(degree: i32) {
    if degree > BSPLIB_MAX_DEGREE as i32 {
        panic!("Standard_OutOfRange: BSplCLib: bspline degree is greater than maximum supported");
    }
}

/// OCCT `BSplCLib_DataContainer_T<1>` (pxx L263-268): the scalar container.
pub(crate) struct DataContainer1 {
    pub poles: [f64; (BSPLIB_MAX_DEGREE + 1) * 2],
    pub knots: [f64; 2 * BSPLIB_MAX_DEGREE],
    pub ders: [f64; 4],
}

impl DataContainer1 {
    pub(crate) fn new() -> Self {
        DataContainer1 {
            poles: [0.0; (BSPLIB_MAX_DEGREE + 1) * 2],
            knots: [0.0; 2 * BSPLIB_MAX_DEGREE],
            ders: [0.0; 4],
        }
    }
}

/// OCCT `BSplCLib_DataContainer_T<2>` (pxx L263-268): the gp_Pnt2d container.
pub(crate) struct DataContainer2 {
    pub poles: [f64; (BSPLIB_MAX_DEGREE + 1) * 3],
    pub knots: [f64; 2 * BSPLIB_MAX_DEGREE],
    pub ders: [f64; 8],
}

impl DataContainer2 {
    pub(crate) fn new() -> Self {
        DataContainer2 {
            poles: [0.0; (BSPLIB_MAX_DEGREE + 1) * 3],
            knots: [0.0; 2 * BSPLIB_MAX_DEGREE],
            ders: [0.0; 8],
        }
    }
}

/// OCCT BSplCLib::LocateParameter(Degree, Knots, Mults, U, Periodic,
/// KnotIndex, NewU) — the 7-argument default-range overload (BSplCLib.cxx
/// L321-364).
pub(crate) fn locate_parameter_span(
    degree: usize,
    knots: &[f64],
    mults: &[i32],
    u: f64,
    periodic: bool,
    knot_index: &mut i32,
    new_u: &mut f64,
) {
    let k_upper = knots.len() as i32;
    let (first, last) = if periodic {
        (1i32, k_upper)
    } else {
        (
            first_uknot_index_mults(degree, mults),
            last_uknot_index_mults(degree, mults),
        )
    };
    if *knot_index < first || *knot_index > last {
        let ufirst = at(knots, first);
        let ulast = at(knots, last);
        locate_parameter_main(knots, u, periodic, first, last, knot_index, new_u, ufirst, ulast);
    } else {
        *new_u = u;
    }
}

/// OCCT BSplCLib::BuildKnots(Degree, Index, Periodic, Knots, Mults, LK)
/// (BSplCLib.cxx L1555-1740).  `Mults` is never null on this path, so only
/// the flat `default:` case (L1690-1701) and the `else` mults branch
/// (L1702-1740) are transcribed; the per-degree unrolled flat cases 1..6
/// (L1570-1689) are byte-for-byte the same loop as the `default` body.
pub(crate) fn build_knots(
    degree: usize,
    index: i32,
    periodic: bool,
    knots: &[f64],
    mults: &[i32],
    knot: &mut [f64],
) {
    let degree_i = degree as i32;
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
    for i in 0..degree as i32 {
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
                knot[(degree_i - 1 - i) as usize] = at(knots, ilow) - loffset;
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
                knot[(degree_i + i) as usize] = at(knots, iupp) + uoffset;
            }
        }
    }
}

/// OCCT BSplCLib::IsRational(Weights, I1, I2, Epsilon = 0) (BSplCLib.cxx
/// L842-861): False if all the weights between I1 and I2 are identical.
pub(crate) fn is_rational_window(weights: &[f64], i1: i32, i2: i32) -> bool {
    let l = weights.len() as i32;
    let i3 = i2 - 1;
    let mut i = i1 - 1;
    while i < i3 {
        if weights[(i % l) as usize] != weights[((i + 1) % l) as usize] {
            return true;
        }
        i += 1;
    }
    false
}

/// OCCT BSplCLib::BuildEval(Degree, Index, Poles, Weights, LP) — the scalar
/// overload (BSplCLib_2.cxx L49-89).  `weights == None` is
/// `BSplCLib::NoWeights()`.
pub(crate) fn build_eval_scalar(
    degree: usize,
    index: i32,
    poles: &[f64],
    weights: Option<&[f64]>,
    lp: &mut [f64],
) {
    let p_upper = poles.len() as i32;
    let mut ip = index; // PLower + Index - 1 with PLower == 1
    match weights {
        None => {
            let mut pole = 0usize;
            for _i in 0..=degree {
                ip += 1;
                if ip > p_upper {
                    ip = 1;
                }
                lp[pole] = at(poles, ip);
                pole += 1;
            }
        }
        Some(w) => {
            let mut pole = 0usize;
            for _i in 0..=degree {
                ip += 1;
                if ip > p_upper {
                    ip = 1;
                }
                let weight = at(w, ip);
                lp[pole + 1] = weight;
                lp[pole] = at(poles, ip) * weight;
                pole += 2;
            }
        }
    }
}

/// OCCT BSplCLib_BuildEval<Point, Vector, Array1OfPoints, 2> (pxx L720-762)
/// — the gp_Pnt2d specialization.
pub(crate) fn build_eval_point2d(
    degree: usize,
    index: i32,
    poles: &[glam::DVec2],
    weights: Option<&[f64]>,
    lp: &mut [f64],
) {
    let p_upper = poles.len() as i32;
    let mut ip = index;
    match weights {
        None => {
            let mut pole = 0usize;
            for _i in 0..=degree {
                ip += 1;
                if ip > p_upper {
                    ip = 1;
                }
                let p = poles[(ip - 1) as usize];
                lp[pole] = p.x;
                lp[pole + 1] = p.y;
                pole += 2;
            }
        }
        Some(w) => {
            let mut pole = 0usize;
            for _i in 0..=degree {
                ip += 1;
                if ip > p_upper {
                    ip = 1;
                }
                let p = poles[(ip - 1) as usize];
                let weight = at(w, ip);
                lp[pole + 2] = weight;
                lp[pole] = p.x * weight;
                lp[pole + 1] = p.y * weight;
                pole += 3;
            }
        }
    }
}

/// OCCT BSplCLib::Bohm — the `case 1:` body (BSplCLib.cxx L1220-1269).
/// `knots` is the 0-based local `knot[2*Degree]` window; `poles` holds
/// `(Degree+1)` scalar rows transformed in place into the value followed by
/// the first `min(N, Degree)` factorial-scaled derivatives.
pub(crate) fn bohm_1d(u: f64, degree: i32, n: i32, knots: &[f64], poles: &mut [f64]) {
    let min = if n < degree { n } else { degree };
    let degm1 = degree - 1;
    let mut ddmi = (degree << 1) + 1;

    // First phase independent of U, compute the poles of the derivatives.
    // psDD = psav + Degree, psDDmDim = psDD - 1 (isize pointer offsets).
    let ps_dd = degree as isize;
    let ps_ddm_dim = ps_dd - 1;
    for i in 0..degree {
        ddmi -= 1;
        let mut pole = ps_dd;
        let mut tbis = ps_ddm_dim;
        let mut jdmi = ddmi;
        let mut j = degm1;
        while j >= i {
            jdmi -= 1;
            poles[pole as usize] -= poles[tbis as usize];
            poles[pole as usize] = if knots[jdmi as usize] == knots[j as usize] {
                0.0
            } else {
                poles[pole as usize] / (knots[jdmi as usize] - knots[j as usize])
            };
            pole -= 1;
            tbis -= 1;
            j -= 1;
        }
    }
    // Second phase, dependant of U.
    let mut idim = -1i32;
    for i in 0..degree {
        idim += 1;
        let mut pole = idim as isize;
        let mut tbis = pole + 1;
        let coef = u - knots[i as usize];
        let mut j = i;
        while j >= 0 {
            poles[pole as usize] += coef * poles[tbis as usize];
            pole -= 1;
            tbis -= 1;
            j -= 1;
        }
    }
    // multiply by the degrees.
    let mut coef = degree as f64;
    let mut dmi = degree;
    let mut pole = 1usize; // psav + 1
    for _i in 1..=min {
        poles[pole] *= coef;
        pole += 1;
        dmi -= 1;
        coef *= dmi as f64;
    }
}

/// OCCT BSplCLib::Bohm — the `case 2:` body (BSplCLib.cxx L1271-1330).
pub(crate) fn bohm_2d(u: f64, degree: i32, n: i32, knots: &[f64], poles: &mut [f64]) {
    let min = if n < degree { n } else { degree };
    let degm1 = degree - 1;
    let mut ddmi = (degree << 1) + 1;

    // First phase independent of U, compute the poles of the derivatives.
    let ps_dd = (degree << 1) as isize;
    let ps_ddm_dim = ps_dd - 2;
    for i in 0..degree {
        ddmi -= 1;
        let mut pole = ps_dd;
        let mut tbis = ps_ddm_dim;
        let mut jdmi = ddmi;
        let mut j = degm1;
        while j >= i {
            jdmi -= 1;
            let coef = if knots[jdmi as usize] == knots[j as usize] {
                0.0
            } else {
                1.0 / (knots[jdmi as usize] - knots[j as usize])
            };
            poles[pole as usize] -= poles[tbis as usize];
            poles[pole as usize] *= coef;
            pole += 1;
            tbis += 1;
            poles[pole as usize] -= poles[tbis as usize];
            poles[pole as usize] *= coef;
            pole -= 3;
            tbis -= 3;
            j -= 1;
        }
    }
    // Second phase, dependant of U.
    let mut idim = -2i32;
    for i in 0..degree {
        idim += 2;
        let mut pole = idim as isize;
        let mut tbis = pole + 2;
        let coef = u - knots[i as usize];
        let mut j = i;
        while j >= 0 {
            poles[pole as usize] += coef * poles[tbis as usize];
            pole += 1;
            tbis += 1;
            poles[pole as usize] += coef * poles[tbis as usize];
            pole -= 3;
            tbis -= 3;
            j -= 1;
        }
    }
    // multiply by the degrees.
    let mut coef = degree as f64;
    let mut dmi = degree;
    let mut pole = 2usize; // psav + 2
    for _i in 1..=min {
        poles[pole] *= coef;
        pole += 1;
        poles[pole] *= coef;
        pole += 1;
        dmi -= 1;
        coef *= dmi as f64;
    }
}

/// OCCT BSplCLib::Bohm — the `case 3:` body (BSplCLib.cxx L1332-1400).
pub(crate) fn bohm_3d(u: f64, degree: i32, n: i32, knots: &[f64], poles: &mut [f64]) {
    let min = if n < degree { n } else { degree };
    let degm1 = degree - 1;
    let mut ddmi = (degree << 1) + 1;

    // First phase independent of U, compute the poles of the derivatives.
    let ps_dd = ((degree << 1) + degree) as isize;
    let ps_ddm_dim = ps_dd - 3;
    for i in 0..degree {
        ddmi -= 1;
        let mut pole = ps_dd;
        let mut tbis = ps_ddm_dim;
        let mut jdmi = ddmi;
        let mut j = degm1;
        while j >= i {
            jdmi -= 1;
            let coef = if knots[jdmi as usize] == knots[j as usize] {
                0.0
            } else {
                1.0 / (knots[jdmi as usize] - knots[j as usize])
            };
            poles[pole as usize] -= poles[tbis as usize];
            poles[pole as usize] *= coef;
            pole += 1;
            tbis += 1;
            poles[pole as usize] -= poles[tbis as usize];
            poles[pole as usize] *= coef;
            pole += 1;
            tbis += 1;
            poles[pole as usize] -= poles[tbis as usize];
            poles[pole as usize] *= coef;
            pole -= 5;
            tbis -= 5;
            j -= 1;
        }
    }
    // Second phase, dependant of U.
    let mut idim = -3i32;
    for i in 0..degree {
        idim += 3;
        let mut pole = idim as isize;
        let mut tbis = pole + 3;
        let coef = u - knots[i as usize];
        let mut j = i;
        while j >= 0 {
            poles[pole as usize] += coef * poles[tbis as usize];
            pole += 1;
            tbis += 1;
            poles[pole as usize] += coef * poles[tbis as usize];
            pole += 1;
            tbis += 1;
            poles[pole as usize] += coef * poles[tbis as usize];
            pole -= 5;
            tbis -= 5;
            j -= 1;
        }
    }
    // multiply by the degrees.
    let mut coef = degree as f64;
    let mut dmi = degree;
    let mut pole = 3usize; // psav + 3
    for _i in 1..=min {
        poles[pole] *= coef;
        pole += 1;
        poles[pole] *= coef;
        pole += 1;
        poles[pole] *= coef;
        pole += 1;
        dmi -= 1;
        coef *= dmi as f64;
    }
}

/// OCCT static PrepareEval (BSplCLib_2.cxx L93-138) — the scalar form.
/// Returns `(dim, rational)`; `u` and `index` are updated in place.
pub(crate) fn prepare_eval_scalar(
    u: &mut f64,
    index: &mut i32,
    degree: usize,
    periodic: bool,
    poles: &[f64],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
    dc: &mut DataContainer1,
) -> (usize, bool) {
    // Set the Index.
    locate_parameter_span(degree, knots, mults, *u, periodic, index, u);

    // make the knots.
    build_knots(degree, *index, periodic, knots, mults, &mut dc.knots);
    *index = pole_index(degree, *index, periodic, mults);

    // check truly rational.
    let mut rational = weights.is_some();
    if let Some(w) = weights {
        // WLower = Weights->Lower() + index, Lower == 1.
        let w_lower = 1 + *index;
        rational = is_rational_window(w, w_lower, w_lower + degree as i32);
    }

    // make the poles.
    if rational {
        build_eval_scalar(degree, *index, poles, weights, &mut dc.poles);
        (2, true)
    } else {
        build_eval_scalar(degree, *index, poles, None, &mut dc.poles);
        (1, false)
    }
}

/// OCCT PrepareEval_T<Point, Vector, Array1OfPoints, 2> (pxx L776-830) — the
/// gp_Pnt2d form.
pub(crate) fn prepare_eval_point2d(
    u: &mut f64,
    index: &mut i32,
    degree: usize,
    periodic: bool,
    poles: &[glam::DVec2],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
    dc: &mut DataContainer2,
) -> (usize, bool) {
    // Set the Index.
    locate_parameter_span(degree, knots, mults, *u, periodic, index, u);

    // make the knots.
    build_knots(degree, *index, periodic, knots, mults, &mut dc.knots);
    *index = pole_index(degree, *index, periodic, mults);

    // check truly rational.
    let mut rational = weights.is_some();
    if let Some(w) = weights {
        let w_lower = 1 + *index;
        rational = is_rational_window(w, w_lower, w_lower + degree as i32);
    }

    // make the poles.
    if rational {
        build_eval_point2d(degree, *index, poles, weights, &mut dc.poles);
        (3, true)
    } else {
        build_eval_point2d(degree, *index, poles, None, &mut dc.poles);
        (2, false)
    }
}

/// OCCT BSplCLib::Bohm dispatch (`switch (Dimension)` at BSplCLib.cxx
/// L1218); the local `knots` window is read 0-based inside Bohm.
fn bohm_dispatch(u: f64, degree: i32, n: i32, knots: &[f64], dimension: i32, poles: &mut [f64]) {
    match dimension {
        1 => bohm_1d(u, degree, n, knots, poles),
        2 => bohm_2d(u, degree, n, knots, poles),
        3 => bohm_3d(u, degree, n, knots, poles),
        // OCCT cases 4 and `default` (L1402-...) coverDimensions > 3; the
        // scalar/point consumers of this module only reach dimensions 1..3.
        _ => panic!("BSplCLib::Bohm: unsupported Dimension {}", dimension),
    }
}

/// OCCT BSplCLib::D0 — the scalar overload (BSplCLib_2.cxx L142-167).
pub(crate) fn bspl_clib_d0_scalar(
    u: f64,
    index: i32,
    degree: usize,
    periodic: bool,
    poles: &[f64],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
    p: &mut f64,
) {
    validate_bspline_degree(degree as i32);
    let mut index = index;
    let mut u = u;
    let mut dc = DataContainer1::new();
    let (dim, rational) = prepare_eval_scalar(
        &mut u, &mut index, degree, periodic, poles, weights, knots, mults, &mut dc,
    );
    bspl_clib_eval_inplace(u, degree as i32, &dc.knots, dim, &mut dc.poles);
    if rational {
        *p = dc.poles[0] / dc.poles[1];
    } else {
        *p = dc.poles[0];
    }
}

/// OCCT BSplCLib::D1 — the scalar overload (BSplCLib_2.cxx L171-197).
pub(crate) fn bspl_clib_d1_scalar(
    u: f64,
    index: i32,
    degree: usize,
    periodic: bool,
    poles: &[f64],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
    p: &mut f64,
    v: &mut f64,
) {
    validate_bspline_degree(degree as i32);
    let mut index = index;
    let mut u = u;
    let mut dc = DataContainer1::new();
    let (dim, rational) = prepare_eval_scalar(
        &mut u, &mut index, degree, periodic, poles, weights, knots, mults, &mut dc,
    );
    bohm_dispatch(u, degree as i32, 1, &dc.knots, dim as i32, &mut dc.poles);
    // result = dc.poles; if (rational) { RationalDerivative(...); result =
    // dc.ders; }  P = result[0]; V = result[1];
    if rational {
        // PLib::RationalDerivative(Degree, 1, 1, *dc.poles, *dc.ders) with
        // the All=true default (PLib.hxx L119-124).
        p_lib_rational_derivative(degree as i32, 1, 1, &dc.poles, &mut dc.ders, true);
        *p = dc.ders[0];
        *v = dc.ders[1];
    } else {
        *p = dc.poles[0];
        *v = dc.poles[1];
    }
}

/// OCCT BSplCLib_D0<Point, Vector, Array1OfPoints, 2> (pxx L843-882) — the
/// gp_Pnt2d D0.
pub(crate) fn bspl_clib_d0_point2d(
    u: f64,
    index: i32,
    degree: usize,
    periodic: bool,
    poles: &[glam::DVec2],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
    p: &mut glam::DVec2,
) {
    validate_bspline_degree(degree as i32);
    let mut index = index;
    let mut u = u;
    let mut dc = DataContainer2::new();
    let (dim, rational) = prepare_eval_point2d(
        &mut u, &mut index, degree, periodic, poles, weights, knots, mults, &mut dc,
    );
    bspl_clib_eval_inplace(u, degree as i32, &dc.knots, dim, &mut dc.poles);

    if rational {
        let w = dc.poles[2]; // dc.poles[Dimension]
        // Traits::CoordsToPointScaled: P = coords / w.
        *p = glam::DVec2::new(dc.poles[0] / w, dc.poles[1] / w);
    } else {
        *p = glam::DVec2::new(dc.poles[0], dc.poles[1]);
    }
}

/// OCCT BSplCLib_D1<Point, Vector, Array1OfPoints, 2> (pxx L896-935) — the
/// gp_Pnt2d D1.
pub(crate) fn bspl_clib_d1_point2d(
    u: f64,
    index: i32,
    degree: usize,
    periodic: bool,
    poles: &[glam::DVec2],
    weights: Option<&[f64]>,
    knots: &[f64],
    mults: &[i32],
    p: &mut glam::DVec2,
    v: &mut glam::DVec2,
) {
    validate_bspline_degree(degree as i32);
    let mut index = index;
    let mut u = u;
    let mut dc = DataContainer2::new();
    let (dim, rational) = prepare_eval_point2d(
        &mut u, &mut index, degree, periodic, poles, weights, knots, mults, &mut dc,
    );
    bohm_dispatch(u, degree as i32, 1, &dc.knots, dim as i32, &mut dc.poles);
    if rational {
        // PLib::RationalDerivative(Degree, 1, Dimension, *dc.poles, *dc.ders)
        // with the All=true default; Dimension == 2.
        p_lib_rational_derivative(degree as i32, 1, 2, &dc.poles, &mut dc.ders, true);
        *p = glam::DVec2::new(dc.ders[0], dc.ders[1]);
        *v = glam::DVec2::new(dc.ders[2], dc.ders[3]);
    } else {
        *p = glam::DVec2::new(dc.poles[0], dc.poles[1]);
        *v = glam::DVec2::new(dc.poles[2], dc.poles[3]);
    }
}


/// OCCT PLib::RationalDerivative(Degree, DerivativeRequest, Dimension, Ders,
/// RDers, All) (PLib.cxx L274-558).  `ders` holds `(Degree+1)` rows of
/// `(Dimension+1)` homogeneous derivative values `[u..., v]` per row (row
/// stride `Dimension + 1` in `ders`, row stride `Dimension` in `rders`).
///
/// OCCT aliases `RationalArray` onto `RDers` when `All == true` and onto a
/// local `derivative_storage` otherwise (then copying the requested row out,
/// L428-437 / L546-556).  Rust cannot re-alias `rders` after borrowing it, so
/// the walk always runs on the local storage buffer (byte-identical to what
/// OCCT writes) and the written span is copied to `rders` when `All == true`;
/// the final content is identical.
pub(crate) fn p_lib_rational_derivative(
    degree: i32,
    derivative_request: i32,
    dimension: i32,
    ders: &[f64],
    rders: &mut [f64],
    all: bool,
) {
    if dimension == 3 {
        // OCCT L326-438 (row stride 4 in `ders`, 3 in the output).
        let de_request1 = derivative_request + 1;
        let min_deg_requ = derivative_request.min(degree);
        let mut binomial_array = vec![1.0f64; de_request1 as usize];
        // OCCT L341-344: DimDeRequ1 = (DeRequest1 << 1) + DeRequest1.
        let mut derivative_storage =
            vec![0.0f64; (((de_request1 << 1) + de_request1) + 1) as usize];

        let mut index = 0usize;
        let mut index1: i32;
        let mut index2 = -6i32;
        let mut other_index = 0usize;
        // OCCT L346: Inverse = 1 / PolesArray[3] (v(0)).
        let inverse = 1.0 / ders[3];

        // OCCT L351-389.
        for ii in 0..=min_deg_requ {
            index2 += 3;
            index1 = index2;
            derivative_storage[index] = ders[other_index];
            index += 1;
            other_index += 1;
            derivative_storage[index] = ders[other_index];
            index += 1;
            other_index += 1;
            derivative_storage[index] = ders[other_index];
            index -= 2;
            other_index += 2;

            let mut jj = ii - 1;
            while jj >= 0 {
                let factor = binomial_array[jj as usize] * ders[(((ii - jj) << 2) + 3) as usize];
                derivative_storage[index] -= factor * derivative_storage[index1 as usize];
                index += 1;
                index1 += 1;
                derivative_storage[index] -= factor * derivative_storage[index1 as usize];
                index += 1;
                index1 += 1;
                derivative_storage[index] -= factor * derivative_storage[index1 as usize];
                index -= 2;
                index1 -= 5;
                jj -= 1;
            }

            for jj in (1..=ii).rev() {
                binomial_array[jj as usize] += binomial_array[(jj - 1) as usize];
            }
            derivative_storage[index] *= inverse;
            index += 1;
            derivative_storage[index] *= inverse;
            index += 1;
            derivative_storage[index] *= inverse;
            index += 1;
        }

        // OCCT L391-426 — rows beyond MinDegRequ.
        for ii in (min_deg_requ + 1)..=derivative_request {
            index2 += 3;
            index1 = index2;
            derivative_storage[index] = 0.0;
            index += 1;
            derivative_storage[index] = 0.0;
            index += 1;
            derivative_storage[index] = 0.0;
            index -= 2;

            let mut jj = ii - 1;
            while jj >= ii - min_deg_requ {
                let factor = binomial_array[jj as usize] * ders[(((ii - jj) << 2) + 3) as usize];
                derivative_storage[index] -= factor * derivative_storage[index1 as usize];
                index += 1;
                index1 += 1;
                derivative_storage[index] -= factor * derivative_storage[index1 as usize];
                index += 1;
                index1 += 1;
                derivative_storage[index] -= factor * derivative_storage[index1 as usize];
                index -= 2;
                index1 -= 5;
                jj -= 1;
            }

            for jj in (1..=ii).rev() {
                binomial_array[jj as usize] += binomial_array[(jj - 1) as usize];
            }
            derivative_storage[index] *= inverse;
            index += 1;
            derivative_storage[index] *= inverse;
            index += 1;
            derivative_storage[index] *= inverse;
            index += 1;
        }

        if all {
            rders[..index].copy_from_slice(&derivative_storage[..index]);
        } else {
            // OCCT L428-437.
            let mut dim_de_requ = ((derivative_request << 1) + derivative_request) as usize;
            rders[0] = derivative_storage[dim_de_requ];
            dim_de_requ += 1;
            rders[1] = derivative_storage[dim_de_requ];
            dim_de_requ += 1;
            rders[2] = derivative_storage[dim_de_requ];
        }
    } else {
        // OCCT L441-557 — the generic Dimension branch (rows of `Dimension`
        // in the output, `Dimension + 1` in `ders`).
        let dimension1 = dimension + 1;
        let dimension2 = dimension << 1;
        let de_request1 = derivative_request + 1;
        let min_deg_requ = derivative_request.min(degree);
        let mut binomial_array = vec![1.0f64; de_request1 as usize];
        // OCCT L458: DimDeRequ1 = Dimension * DeRequest1.
        let mut derivative_storage = vec![0.0f64; (dimension * de_request1) as usize];

        let mut index = 0usize;
        let mut index1: i32;
        let mut index2 = -dimension2;
        let mut other_index = 0usize;
        // OCCT L463: Inverse = 1 / PolesArray[Dimension] (v(0)).
        let inverse = 1.0 / ders[dimension as usize];

        // OCCT L468-506.
        for ii in 0..=min_deg_requ {
            index2 += dimension;
            index1 = index2;

            for _kk in 0..dimension {
                derivative_storage[index] = ders[other_index];
                index += 1;
                other_index += 1;
            }
            index -= dimension as usize;
            other_index += 1;

            let mut jj = ii - 1;
            while jj >= 0 {
                let factor =
                    binomial_array[jj as usize] * ders[((ii - jj) * dimension1 + dimension) as usize];
                for _kk in 0..dimension {
                    derivative_storage[index] -= factor * derivative_storage[index1 as usize];
                    index += 1;
                    index1 += 1;
                }
                index -= dimension as usize;
                index1 -= dimension2;
                jj -= 1;
            }

            for jj in (1..=ii).rev() {
                binomial_array[jj as usize] += binomial_array[(jj - 1) as usize];
            }

            for _kk in 0..dimension {
                derivative_storage[index] *= inverse;
                index += 1;
            }
        }

        // OCCT L508-544 — rows beyond MinDegRequ.
        for ii in (min_deg_requ + 1)..=derivative_request {
            index2 += dimension;
            index1 = index2;

            for _kk in 0..dimension {
                derivative_storage[index] = 0.0;
                index += 1;
            }
            index -= dimension as usize;

            let mut jj = ii - 1;
            while jj >= ii - min_deg_requ {
                let factor =
                    binomial_array[jj as usize] * ders[((ii - jj) * dimension1 + dimension) as usize];
                for _kk in 0..dimension {
                    derivative_storage[index] -= factor * derivative_storage[index1 as usize];
                    index += 1;
                    index1 += 1;
                }
                index -= dimension as usize;
                index1 -= dimension2;
                jj -= 1;
            }

            for jj in (1..=ii).rev() {
                binomial_array[jj as usize] += binomial_array[(jj - 1) as usize];
            }

            for _kk in 0..dimension {
                derivative_storage[index] *= inverse;
                index += 1;
            }
        }

        if all {
            rders[..index].copy_from_slice(&derivative_storage[..index]);
        } else {
            // OCCT L546-556.
            let mut dim_de_requ = (dimension * derivative_request) as usize;
            for kk in 0..dimension as usize {
                rders[kk] = derivative_storage[dim_de_requ];
                dim_de_requ += 1;
            }
        }
    }
}
