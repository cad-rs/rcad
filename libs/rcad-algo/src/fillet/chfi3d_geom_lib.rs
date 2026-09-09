//! OCCT TKGeomBase translations consumed by the ElSpine construction
//! (`chfi3d_perform_elspine`): GeomLib::ExtendCurveToPoint /
//! AdjustExtremity and their PLib/math dependencies.
//!
//! Sources (1:1, per-function line anchors below):
//!   - TKGeomBase/GeomLib/GeomLib.cxx L1269-1413 (ExtendCurveToPoint),
//!     L1169-1265 (AdjustExtremity), L126-253 (static ComputeLambda)
//!   - TKMath/PLib/PLib.cxx L1404-1478 (HermiteCoefficients),
//!     L1482-1605 (CoefficientsPoles), PLib::Bin
//!   - TKGeomBase/GeomLib/GeomLib_PolyFunc.cxx L18-60
//!
//! Architecture mappings: math_Matrix -> kernel Matrix, math_Vector ->
//! kernel Vector, math_Gauss -> kernel MathGauss, math::GaussPoints/
//! GaussWeights -> kernel gauss_points/gauss_weights, PLib::EvalPolynomial
//! -> kernel eval_polynomial_flat, PLib::NoDerivativeEvalPolynomial ->
//! kernel no_derivative_eval_polynomial_flat.

use glam::DVec3;
use rcad_kernel::base::convert::ConvertParameterisation;
use rcad_kernel::geom::{BSplineCurve3, Curve3};
use rcad_kernel::math::gauss_points::{gauss_points, gauss_weights};
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::math_matrix::Matrix;
use rcad_kernel::math::math_matrix::Vector;
use rcad_kernel::math::plib::{eval_polynomial_flat, no_derivative_eval_polynomial_flat};
use rcad_kernel::math::root::function_all_roots::{
    FunctionAllRoots, FunctionSample, FunctionValue, FunctionWithDerivative,
};
use rcad_kernel::math::VecD;

use super::chfi3d_perform_elspine::{bspline_value, flat_knots, GeomConvertCompCurveToBSplineCurve};

// =========================================================================
// OCCT FoundationClasses/TKMath/PLib/PLib.cxx — PLib::Bin(N, P) (the binomial
// coefficient as a double).
// =========================================================================
fn plib_bin(n: i32, p: i32) -> f64 {
    let mut res = 1.0f64;
    for i in 1..=p {
        res *= (n - p + i) as f64 / i as f64;
    }
    res
}

// =========================================================================
// OCCT FoundationClasses/TKMath/PLib/PLib.cxx L1482-1493 (gp_Pnt overload) +
// L1522-1605 (generic, dim = 3, non-rational) — PLib::CoefficientsPoles:
// power-basis coefficients -> Bernstein poles of the same degree.
// =========================================================================
pub(crate) fn plib_coefficients_poles_dim3(coefs: &[DVec3]) -> Vec<DVec3> {
    let reflen = coefs.len();
    let mut poles = vec![DVec3::ZERO; reflen];
    // Les Extremites (OCCT L1548-1551).
    poles[0] = coefs[0];
    poles[reflen - 1] = coefs[reflen - 1];
    // OCCT L1566-1576: Poles(i) = Coefs(i) / Bin(reflen-1, i-1).
    for i in 2..reflen as i32 {
        let cnp = plib_bin(reflen as i32 - 1, i - 1);
        poles[(i - 1) as usize] = coefs[(i - 1) as usize] / cnp;
    }
    // OCCT L1579-1591: the Pascal accumulation.
    for i in 1..=reflen as i32 - 1 {
        let mut j = reflen as i32 - 1;
        while j >= i {
            let prev = poles[(j - 1) as usize];
            poles[j as usize] += prev;
            j -= 1;
        }
    }
    poles
}

// =========================================================================
// OCCT FoundationClasses/TKMath/PLib/PLib.cxx L1404-1478 —
// PLib::HermiteCoefficients(FirstParameter, LastParameter, FirstOrder,
// LastOrder, MatrixCoefs): the Hermite coefficient matrix solved by
// math_Gauss.  math_Vector -> Vector, math_Matrix -> Matrix (kernel
// translations); the OCCT two-argument Solve(B, Coeff) maps onto the
// in-place kernel solve of a copy of B.
// =========================================================================
pub(crate) fn plib_hermite_coefficients(
    first_parameter: f64,
    last_parameter: f64,
    first_order: i32,
    last_order: i32,
    matrix_coefs: &mut Matrix,
) -> bool {
    let nb_coeff = first_order + last_order + 2;
    let mut ordre = [0i32; 2];
    let mut iof = 0i32;
    let mut t_borne = first_parameter;
    let mut coeff = Vector::new_init(1, nb_coeff, 0.0);
    let mut b = Vector::new_init(1, nb_coeff, 0.0);
    let mut mat = Matrix::new_init(1, nb_coeff, 1, nb_coeff, 0.0);

    // Test de validites (OCCT L1416-1435).
    if first_order < 0 || last_order < 0 {
        return false;
    }
    let d1 = first_parameter.abs();
    let d2p = last_parameter.abs();
    if d1 > 100.0 || d2p > 100.0 {
        return false;
    }
    let d2 = d2p + d1;
    if d2 < 0.01 {
        return false;
    }
    if (last_parameter - first_parameter).abs() / d2 < 0.01 {
        return false;
    }

    // Calcul de la matrice a inverser (OCCT L1437-1461).
    ordre[0] = first_order + 1;
    ordre[1] = last_order + 1;
    for cote in 0..=1 {
        coeff = Vector::new_init(1, nb_coeff, 1.0);
        for pp in 1..=ordre[cote] {
            let ii = pp + iof;
            let mut prod = 1.0f64;
            for jj in pp..=nb_coeff {
                mat.set(ii, jj, coeff.get(jj) * prod);
                coeff.set(jj, coeff.get(jj) * (jj - pp) as f64);
                prod *= t_borne;
            }
        }
        t_borne = last_parameter;
        iof = ordre[0];
    }

    // resolution du systemes (OCCT L1463-1477).
    let resol_coeff = MathGauss::with_min_pivot(&mat.data, 1.0e-10);
    if !resol_coeff.is_done() {
        return false;
    }
    for ii in 1..=nb_coeff {
        b.set(ii, 1.0);
        // OCCT: ResolCoeff.Solve(B, Coeff) — the kernel solves in place on a
        // copy of B.
        let mut x = b.data.clone();
        resol_coeff.solve(&mut x);
        for jj in 1..=nb_coeff {
            coeff.set(jj, x.get(jj as usize));
        }
        for jj in 1..=nb_coeff {
            matrix_coefs.set(ii, jj, coeff.get(jj));
        }
        b.set(ii, 0.0);
    }
    true
}

// =========================================================================
// OCCT math_Vector operators (math_Vector.lxx) used by ComputeLambda —
// Norm2, the scalar product Multiplied(Vector), and the scaled divide.
// =========================================================================
fn vector_norm2(v: &Vector) -> f64 {
    let mut s = 0.0f64;
    for i in v.lower()..=v.upper() {
        let x = v.get(i);
        s += x * x;
    }
    s
}

fn vector_multiplied(a: &Vector, b: &Vector) -> f64 {
    let mut s = 0.0f64;
    for i in a.lower()..=a.upper() {
        s += a.get(i) * b.get(i);
    }
    s
}

fn vector_divided_scalar(v: &Vector, s: f64) -> Vector {
    let mut r = v.clone();
    for i in r.lower()..=r.upper() {
        r.set(i, r.get(i) / s);
    }
    r
}

// =========================================================================
// OCCT TKGeomBase/GeomLib/GeomLib_PolyFunc.cxx L18-60 — the derivative
// polynomial of the integrated square-deviation polynomial, evaluated
// through PLib::EvalPolynomial.  Implements the kernel
// math_FunctionWithDerivative trait (OCCT math_FunctionWithDerivative).
// =========================================================================
struct GeomLibPolyFunc {
    /// OCCT: math_Vector myCoeffs(1, Coeffs.Length() - 1).
    my_coeffs: Vec<f64>,
}

impl GeomLibPolyFunc {
    /// OCCT GeomLib_PolyFunc::GeomLib_PolyFunc(const math_Vector& Coeffs):
    /// myCoeffs(ii) = ii * Coeffs(ii + 1).
    fn new(coeffs: &[f64]) -> Self {
        let n = coeffs.len() - 1;
        let mut my_coeffs = vec![0.0f64; n];
        for ii in 1..=n {
            my_coeffs[ii - 1] = ii as f64 * coeffs[ii];
        }
        GeomLibPolyFunc { my_coeffs }
    }
}

impl FunctionValue for GeomLibPolyFunc {
    /// OCCT GeomLib_PolyFunc::Value(X, F): EvalPolynomial(X, 0, Len-1, 1).
    fn value(&mut self, x: f64) -> Option<f64> {
        let mut f = [0.0f64; 1];
        eval_polynomial_flat(
            x,
            0,
            self.my_coeffs.len() as i32 - 1,
            1,
            &self.my_coeffs,
            &mut f,
        );
        Some(f[0])
    }
}

impl FunctionWithDerivative for GeomLibPolyFunc {
    /// OCCT GeomLib_PolyFunc::Derivative(X, D): EvalPolynomial(X, 1, Len-1, 1).
    fn derivative(&mut self, x: f64) -> Option<f64> {
        let mut aux = [0.0f64; 2];
        eval_polynomial_flat(
            x,
            1,
            self.my_coeffs.len() as i32 - 1,
            1,
            &self.my_coeffs,
            &mut aux,
        );
        Some(aux[1])
    }

    /// OCCT GeomLib_PolyFunc::Values(X, F, D).
    fn values(&mut self, x: f64) -> Option<(f64, f64)> {
        let mut aux = [0.0f64; 2];
        eval_polynomial_flat(
            x,
            1,
            self.my_coeffs.len() as i32 - 1,
            1,
            &self.my_coeffs,
            &mut aux,
        );
        Some((aux[0], aux[1]))
    }
}

// =========================================================================
// OCCT TKGeomBase/GeomLib/GeomLib.cxx L126-253 — the static ComputeLambda
// (Constraint, Hermit, Length, Lambda): the Gauss-quadrature minimization of
// the integrated square deviation |C'(t)|^2 - |C'(0)|^2 over [0, 1],
// followed by the extremum refinement.
// GAP: the extremum refinement runs OCCT math_FunctionAllRoots over a
// GeomLib_LogSample (log-spaced parameters); the kernel math_FunctionAllRoots
// binds the concrete linear math_FunctionSample (no virtual GetParameter),
// so the refinement samples the same interval linearly — pending the TKMath
// batch.  The OCCT fallback (keep Lambda when no better extremum is found)
// is preserved.
// =========================================================================
fn geom_lib_compute_lambda(constraint: &Matrix, hermit: &Matrix, length: f64, lambda: &mut f64) {
    let size = hermit.row_number() as usize;
    let continuity = size - 2;

    // Minimization (OCCT L135-143): HDer(ii, jj) = ii * Hermit(jj, ii + 1),
    // flattened row-major (the &HDer(1, 1) view).
    let mut hder = vec![0.0f64; (size - 1) * size];
    for jj in 1..=size as i32 {
        for ii in 1..size as i32 {
            hder[((ii - 1) * size as i32 + (jj - 1)) as usize] =
                ii as f64 * hermit.get(jj, ii + 1);
        }
    }

    let r = constraint.row_number() as usize;
    let mut vec1 = constraint.col(1);
    let mut vec2 = constraint.col(2);
    // OCCT L148-149: Vec3/Vec4 are declared with the row size only and stay
    // untouched for Continuity < 2 / < 3.
    let mut vec3 = Vector::new(1, r as i32);
    let mut vec4 = Vector::new(1, r as i32);
    let mut v = vec![0.0f64; size];

    // OCCT L153-155: Vec2 = Constraint.Col(2); Vec2 /= Length; squared1.
    vec2 = vector_divided_scalar(&vec2, length);
    let squared1 = vector_norm2(&vec2);

    // OCCT L161-165: the Gauss quadrature.
    let g_ordre = 4 + 4 * continuity as i32;
    let d_dim = (continuity * (continuity + 2)) as i32;
    let mut gauss_p = VecD::new(g_ordre as usize);
    let mut gauss_w = VecD::new(g_ordre as usize);
    gauss_points(g_ordre as usize, &mut gauss_p);
    gauss_weights(g_ordre as usize, &mut gauss_w);
    let mut pol2 = vec![0.0f64; (2 * continuity + 1) as usize];
    let mut pol4 = vec![0.0f64; (4 * continuity + 1) as usize];

    for ip in 1..=g_ordre {
        let t = (gauss_p.get(ip as usize) + 1.0) / 2.0;
        let gw = gauss_w.get(ip as usize);
        // OCCT L172: PLib::NoDerivativeEvalPolynomial(t, Continuity,
        // Continuity + 2, DDim, polynome[0], valhder[0]).
        no_derivative_eval_polynomial_flat(
            t,
            continuity as i32,
            size as i32,
            d_dim,
            &hder,
            &mut v,
        );
        for x in v.iter_mut() {
            *x /= length; // OCCT L173: V /= Length (Normalisation)
        }

        //                      i
        // C'(t) = SUM Vi*Lambda  (OCCT L175-191).
        for i in 1..=r as i32 {
            vec1.set(
                i,
                constraint.get(i, 1) * v[0] + v[size - 1] * constraint.get(i, size as i32),
            );
        }
        for i in 1..=r as i32 {
            vec2.set(i, constraint.get(i, 2) * v[1]);
        }
        if continuity > 1 {
            for i in 1..=r as i32 {
                vec3.set(i, constraint.get(i, 3) * v[2]);
            }
            if continuity > 2 {
                for i in 1..=r as i32 {
                    vec4.set(i, constraint.get(i, 4) * v[3]);
                }
            }
        }

        //   2          2
        // C'(t) - C'(0)   (OCCT L193-211).
        pol2[0] = vector_norm2(&vec1);
        pol2[1] = 2.0 * vector_multiplied(&vec1, &vec2);
        pol2[2] = vector_norm2(&vec2) - squared1;
        if continuity > 1 {
            pol2[2] += 2.0 * vector_multiplied(&vec1, &vec3);
            pol2[3] = 2.0 * vector_multiplied(&vec2, &vec3);
            pol2[4] = vector_norm2(&vec3);
            if continuity > 2 {
                pol2[3] += 2.0 * vector_multiplied(&vec1, &vec4);
                pol2[4] += 2.0 * vector_multiplied(&vec2, &vec4);
                pol2[5] = 2.0 * vector_multiplied(&vec3, &vec4);
                pol2[6] = vector_norm2(&vec4);
            }
        }

        //                     2      2  2
        // Integrale de ( C'(t) - C'(0) )  (OCCT L213-223).
        for ii in 1..=pol2.len() {
            let mut pp = ii;
            for jj in 1..ii {
                pol4[pp - 1] += 2.0 * gw * pol2[ii - 1] * pol2[jj - 1];
                pp += 1;
            }
            pol4[2 * ii - 2] += gw * pol2[ii - 1] * pol2[ii - 1];
        }
    }

    let mut e_min = 0.0f64;
    no_derivative_eval_polynomial_flat(
        *lambda,
        pol4.len() as i32 - 1,
        1,
        pol4.len() as i32 - 1,
        &pol4,
        std::slice::from_mut(&mut e_min),
    );

    if e_min > rcad_kernel::core::precision::CONFUSION {
        // Search for extrema of the function (OCCT L229-252).
        let mut ff = GeomLibPolyFunc::new(&pol4);
        let sample = FunctionSample::new(*lambda / 1000.0, 50.0 * *lambda, 100);
        let solve = FunctionAllRoots::new(
            &mut ff,
            &sample,
            rcad_kernel::core::precision::CONFUSION,
            rcad_kernel::core::precision::CONFUSION * (length + 1.0),
            1.0e-15,
        );
        if solve.is_done() {
            for ii in 1..=solve.nb_points() {
                let t = solve.get_point(ii);
                let mut e = 0.0f64;
                no_derivative_eval_polynomial_flat(
                    t,
                    pol4.len() as i32 - 1,
                    1,
                    pol4.len() as i32 - 1,
                    &pol4,
                    std::slice::from_mut(&mut e),
                );
                if e < e_min {
                    *lambda = t;
                    e_min = e;
                }
            }
        }
    }
}

// =========================================================================
// OCCT TKGeomBase/GeomLib/GeomLib.cxx L1269-1413 — ExtendCurveToPoint
// (Curve, Point, Continuity, After): extend a bounded BSpline to `point`
// with the given continuity at the join, through a degree-(Continuity+2)
// Bezier extension concatenated by GeomConvert_CompCurveToBSplineCurve
// (knot-vector consistent by construction).  (This replaces the former
// kernel extend_curve_to_point carrier, whose Start/End branch dropped
// `n_first_min - 1` clamped knot repeats while adding a single pole,
// breaking sum(mults) = nbPoles + degree + 1.)
// =========================================================================
pub(crate) fn geom_lib_extend_curve_to_point(
    curve: &mut BSplineCurve3,
    point: DVec3,
    continuity: i32,
    after: bool,
) {
    if continuity < 1 || continuity > 3 {
        return;
    }
    let size = (continuity + 2) as usize;
    let mut tol = 1.0e-6;
    let mut lambda;
    // Convert the input (OCCT L1286): GeomConvert_CompCurveToBSplineCurve
    // Concat(Curve, Convert_QuasiAngular) — the current curve is already a
    // BSpline, the ctor copies it.
    let mut concat = GeomConvertCompCurveToBSplineCurve::with_basis_curve(
        &Curve3::BSpline(curve.clone()),
        ConvertParameterisation::QuasiAngular,
    );

    // Construction constraints (OCCT L1288-1309).
    let ubord = if after {
        curve.last_parameter()
    } else {
        curve.first_parameter()
    };
    let mut mat_coefs = Matrix::new_init(1, size as i32, 1, size as i32, 0.0);
    // OCCT ignores the HermiteCoefficients return value here (L1298).
    let _ = plib_hermite_coefficients(0.0, 1.0, continuity, 0, &mut mat_coefs);

    // OCCT L1304: Curve->D3(Ubord, p0, d1, d2, d3).
    let p0 = bspline_value(curve, ubord);
    let mut d1 = curve.dn(ubord, 1);
    let d2 = curve.dn(ubord, 2);
    let mut d3 = curve.dn(ubord, 3);
    if !after {
        // Invert the parameterization (OCCT L1305-1309: d1 and d3 only).
        d1 = -d1;
        d3 = -d3;
    }

    let l1 = p0.distance(point);
    if l1 > tol {
        // Lambda is the ratio to apply to the derivative of the curve to
        // obtain the derivative of the extension (OCCT L1311-1343).
        let f = curve.first_parameter();
        let mut dt = (curve.last_parameter() - f) / 9.0;
        let mut norm = d1.length();
        let mut t = f + dt;
        for _ii in 1..=8i32 {
            let daux = curve.dn(t, 1);
            norm += daux.length();
            t += dt;
        }
        norm /= 9.0;
        dt = d1.length() / norm;
        if dt < 1.5 && dt > 0.75 {
            // The edge is within the average, keep it
            lambda = 1.0 / (d1.length() / l1).max(tol);
        } else {
            lambda = 1.0 / (norm / l1).max(tol);
        }
    } else {
        return; // No extension
    }

    // Optimization of Lambda (OCCT L1345-1368).
    let mut cons = Matrix::new_init(1, 3, 1, size as i32, 0.0);
    cons.set(1, 1, p0.x);
    cons.set(2, 1, p0.y);
    cons.set(3, 1, p0.z);
    cons.set(1, 2, d1.x);
    cons.set(2, 2, d1.y);
    cons.set(3, 2, d1.z);
    cons.set(1, size as i32, point.x);
    cons.set(2, size as i32, point.y);
    cons.set(3, size as i32, point.z);
    if continuity >= 2 {
        cons.set(1, 3, d2.x);
        cons.set(2, 3, d2.y);
        cons.set(3, 3, d2.z);
    }
    if continuity >= 3 {
        cons.set(1, 4, d3.x);
        cons.set(2, 4, d3.y);
        cons.set(3, 4, d3.z);
    }
    geom_lib_compute_lambda(&cons, &mat_coefs, l1, &mut lambda);

    // Construction in the Polynomial Basis (OCCT L1370-1394).
    let mut cont = vec![DVec3::ZERO; size];
    cont[0] = p0;
    cont[1] = d1 * lambda;
    if continuity >= 2 {
        cont[2] = d2 * lambda.powi(2);
    }
    if continuity >= 3 {
        cont[3] = d3 * lambda.powi(3);
    }
    cont[size - 1] = point;

    let mut extra_coeffs = vec![DVec3::ZERO; size];
    for ii in 1..=size as i32 {
        for jj in 1..=size as i32 {
            extra_coeffs[(jj - 1) as usize] += mat_coefs.get(ii, jj) * cont[(ii - 1) as usize];
        }
    }

    // Conversion to the Bernstein Basis (OCCT L1396-1397).
    let extra_poles = plib_coefficients_poles_dim3(&extra_coeffs);

    // OCCT L1399: Geom_BezierCurve(ExtrapPoles) — a Bezier is the
    // single-span BSpline with knots [0, 1], mults [size, size].
    let bezier = Curve3::BSpline(BSplineCurve3 {
        degree: size - 1,
        knots: flat_knots(&[0.0, 1.0], &[size as i32, size as i32]),
        control_points: extra_poles.clone(),
        weights: vec![1.0; size],
        is_periodic: false,
    });

    let dist = extra_poles[0].distance(p0);
    tol += dist;

    // Concatenation (OCCT L1405-1410).  The OCCT Add downcasts the
    // Geom_BezierCurve to Geom_BSplineCurve (null) and reconverts through
    // CurveToBSplineCurve, which for a Bezier returns the identical
    // single-span BSpline — the direct downcast-copy path of Add is
    // algebraically identical.
    let ok = concat.add(&bezier, tol, after);
    if !ok {
        panic!("Standard_ConstructionError: ExtendCurveToPoint");
    }

    *curve = concat
        .bspline_curve()
        .expect("ExtendCurveToPoint: null concatenated curve");
}

// =========================================================================
// OCCT TKGeomBase/GeomLib/GeomLib.cxx L1169-1265 — AdjustExtremity.
// GAP: the Hermite pole deformation needs PLib::CoefficientsPoles on
// rational pole arrays and Geom_BSplineCurve::InsertKnot — the closest
// translated kernel machinery is Geom_BSplineCurve::MovePointAndTangent
// (kernel BSplineCurve3, the OCCT L1123-1158 translation), applied at both
// extremities.  Pole-count preserving, so the knot-vector invariant holds.
// =========================================================================
pub(crate) fn geom_lib_adjust_extremity(
    curve: &mut BSplineCurve3,
    p1: DVec3,
    p2: DVec3,
    t1: DVec3,
    t2: DVec3,
) {
    let first = curve.first_parameter();
    let last = curve.last_parameter();
    let _ = curve.move_point_and_tangent(first, p1, t1.normalize_or_zero(), 1.0e-7, -1, 1);
    let _ = curve.move_point_and_tangent(last, p2, t2.normalize_or_zero(), 1.0e-7, -1, 1);
}
