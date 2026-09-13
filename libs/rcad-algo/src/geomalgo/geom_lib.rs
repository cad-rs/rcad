//! OCCT GeomLib statics (TKGeomBase/GeomLib) — the 1:1 translation of
//! `GeomLib::ExtendCurveToPoint` (GeomLib.cxx L1269-1410) and of the file
//! static `ComputeLambda` (GeomLib.cxx L126-253) that it calls, of
//! `GeomLib::To3d` (GeomLib.cxx L559-675) and of the curve-on-surface
//! approximation entry point `GeomLib::BuildCurve3d` (GeomLib.cxx L1051-1163)
//! with its evaluator `GeomLib_CurveOnSurfaceEvaluator` (GeomLib.cxx
//! L974-1050).
//!
//! `ExtendCurveToPoint(Curve, Point, Continuity, After)` extends the bounded
//! curve to `Point` with a degree-(Continuity+2) Bezier segment built from
//! the Hermite constraints at the join and concatenated by
//! `GeomConvert_CompCurveToBSplineCurve` (which preserves the join to the
//! requested tolerance).  The rcad signature takes `&mut Curve3` for the OCCT
//! `Handle(Geom_BoundedCurve)&` — the rcad `Curve3` variants
//! BSpline / Bezier / Trimmed are the OCCT `Geom_BoundedCurve`
//! instantiations, and `CurveEval::{point_at, derivative_at, derivative2_at,
//! derivative3_at}` are the OCCT `Geom_Curve::D3(U, P, D1, D2, D3)` virtual
//! dispatch.
//!
//! `BuildCurve3d` takes rcad's kernel `CurveOnSurface` (the real
//! `Adaptor3d_CurveOnSurface`) as the input adaptor; the two OCCT down-casts
//! to `GeomAdaptor_Surface` / `Geom2dAdaptor_Curve` ride the
//! `kernel_surface` / `kernel_curve2d` bridges, and the null-handle outcome of
//! a failed down-cast is their `None`.  The isoline branch
//! (`isIsoLine` / `buildC3dOnIsoLine`) and the result builder
//! (`GeomLib_MakeCurvefromApprox`) are the sibling translations
//! `geomalgo::geom_lib_iso_line` / `geomalgo::geom_lib_make_curve_from_approx`.
//!
//! The re-used translations (outside this module's OCCT class):
//! - `GeomConvert_CompCurveToBSplineCurve` (TKGeomBase/GeomConvert):
//!   `crate::fillet::chfi3d_perform_elspine::GeomConvertCompCurveToBSplineCurve`;
//! - `PLib::HermiteCoefficients` / `PLib::CoefficientsPoles`:
//!   `crate::fillet::chfi3d_geom_lib::{plib_hermite_coefficients,
//!   plib_coefficients_poles_dim3}`;
//! - `math::GaussPoints` / `math::GaussWeights` /
//!   `PLib::NoDerivativeEvalPolynomial` / `math_FunctionAllRoots`:
//!   the rcad-kernel math translations.
//!
//! Known GAP (shared with the `chfi3d_geom_lib.rs` translation of the same
//! OCCT static): the extremum refinement of `ComputeLambda` runs
//! `math_FunctionAllRoots` over a `GeomLib_LogSample` (log-spaced
//! parameters); the kernel `math_FunctionAllRoots` binds the linear
//! `math_FunctionSample`, so the refinement samples the same interval
//! linearly — the OCCT fallback (keep Lambda when no better extremum is
//! found) is preserved.

use std::sync::Arc;

use glam::DVec3;
use rcad_kernel::base::convert::ConvertParameterisation;
use rcad_kernel::base::proj_lib::adaptor::{Adaptor3dCurve, Curve2dHandle, CurveOnSurface};
use rcad_kernel::core::precision::p_confusion;
use rcad_kernel::geom::{BezierCurve3, Curve2d, Curve3, CurveEval, Surface3};
use rcad_kernel::math::adv_approx::{ApproxAFunction, EvaluatorFunction, PrefAndRec};
use rcad_kernel::math::gauss_points::{gauss_points, gauss_weights};
use rcad_kernel::math::math_matrix::{Matrix, Vector};
use rcad_kernel::math::plib::{eval_polynomial_flat, no_derivative_eval_polynomial_flat};
use rcad_kernel::math::root::function_all_roots::{
    FunctionAllRoots, FunctionSample, FunctionValue, FunctionWithDerivative,
};
use rcad_kernel::math::GeomAbsShape;
use rcad_kernel::math::VecD;

use crate::fillet::chfi3d_geom_lib::{plib_coefficients_poles_dim3, plib_hermite_coefficients};
use crate::fillet::chfi3d_perform_elspine::GeomConvertCompCurveToBSplineCurve;
use crate::geomalgo::geom_lib_iso_line::{build_c3d_on_iso_line, is_iso_line};
use crate::geomalgo::geom_lib_make_curve_from_approx::GeomLibMakeCurvefromApprox;

/// OCCT GeomLib_PolyFunc (GeomLib.hxx L54-90) — the polynomial whose roots
/// the extremum refinement of ComputeLambda searches: the rcad
/// `math_FunctionWithDerivative` binding of the OCCT class.
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

/// OCCT static ComputeLambda(Constraint, Hermit, Length, Lambda)
/// (GeomLib.cxx L126-253) — the Gauss-quadrature minimization of the
/// integrated square deviation `|C'(t)|^2 - |C'(0)|^2` over [0, 1] followed
/// by the extremum refinement.  See the module header for the LogSample GAP.
fn compute_lambda(constraint: &Matrix, hermit: &Matrix, length: f64, lambda: &mut f64) {
    let size = hermit.row_number() as usize;
    let continuity = size - 2;

    // Minimization (OCCT L135-143): HDer(ii, jj) = ii * Hermit(jj, ii + 1) —
    // the flattened &HDer(1, 1) view the polynomial evaluator reads.
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
    // OCCT L148-149: Vec3/Vec4 are declared with the row size only, and stay
    // untouched for Continuity < 2 / < 3.
    let mut vec3 = Vector::new_init(1, r as i32, 0.0);
    let mut vec4 = Vector::new_init(1, r as i32, 0.0);
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
        no_derivative_eval_polynomial_flat(t, continuity as i32, size as i32, d_dim, &hder, &mut v);
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

    // OCCT L225-226: PLib::NoDerivativeEvalPolynomial(Lambda,
    // pol4.Length() - 1, 1, pol4.Length() - 1, pol4(1), EMin).
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
        // Search for extrema of the function (OCCT L229-252): the
        // GeomLib_PolyFunc over the GeomLib_LogSample window, refined by
        // math_FunctionAllRoots.
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

/// OCCT math_Vector::Norm2() (math_Vector.hxx).
fn vector_norm2(v: &Vector) -> f64 {
    let mut s = 0.0f64;
    for i in 1..=v.length() {
        s += v.get(i) * v.get(i);
    }
    s
}

/// OCCT math_Vector::Multiplied(Vector) (math_Vector.cxx).
fn vector_multiplied(a: &Vector, b: &Vector) -> f64 {
    let mut s = 0.0f64;
    for i in 1..=a.length() {
        s += a.get(i) * b.get(i);
    }
    s
}

/// OCCT math_Vector::operator/=(double) (math_Vector.cxx).
fn vector_divided_scalar(v: &Vector, s: f64) -> Vector {
    let mut out = Vector::new_init(1, v.length(), 0.0);
    for i in 1..=v.length() {
        out.set(i, v.get(i) / s);
    }
    out
}

/// OCCT Geom_Curve::D3(U, P, D1, D2, D3) over the rcad `Curve3` — the
/// virtual dispatch is the `CurveEval` one-derivative-per-call surface.
fn curve_d3(the_curve: &Curve3, the_u: f64) -> (DVec3, DVec3, DVec3, DVec3) {
    (
        CurveEval::point_at(the_curve, the_u),
        CurveEval::derivative_at(the_curve, the_u),
        CurveEval::derivative2_at(the_curve, the_u),
        CurveEval::derivative3_at(the_curve, the_u),
    )
}

/// OCCT GeomLib::ExtendCurveToPoint(Curve, Point, Continuity, After)
/// (GeomLib.cxx L1269-1410).
pub fn extend_curve_to_point(
    the_curve: &mut Curve3,
    the_point: DVec3,
    the_continuity: i32,
    the_after: bool,
) {
    // OCCT L1276-1279.
    if the_continuity < 1 || the_continuity > 3 {
        return;
    }
    let size = the_continuity + 2;
    let mut tol = 1.0e-6;
    let mut lambda = 0.0f64;

    // OCCT L1286: GeomConvert_CompCurveToBSplineCurve Concat(Curve,
    // Convert_QuasiAngular).
    let mut concat = GeomConvertCompCurveToBSplineCurve::with_basis_curve(
        the_curve,
        ConvertParameterisation::QuasiAngular,
    );

    // Construction constraints (OCCT L1288-1298).
    let a_dom = CurveEval::default_domain(the_curve);
    let ubord = if the_after { a_dom[1] } else { a_dom[0] };
    let mut mat_coefs = Matrix::new_init(1, size as i32, 1, size as i32, 0.0);
    // OCCT L1298 ignores the PLib::HermiteCoefficients return value.
    let _ = plib_hermite_coefficients(0.0, 1.0, the_continuity, 0, &mut mat_coefs);

    // OCCT L1304: Curve->D3(Ubord, p0, d1, d2, d3).
    let (p0, mut d1, d2, mut d3) = curve_d3(the_curve, ubord);
    if !the_after {
        // OCCT L1305-1309: the parameterization inversion (d1 and d3 only).
        d1 = -d1;
        d3 = -d3;
    }

    // OCCT L1311-1343: L1 and the Lambda ratio.
    let l1 = p0.distance(the_point);
    if l1 > tol {
        let f = a_dom[0];
        let mut dt = (a_dom[1] - f) / 9.0;
        let mut norm = d1.length();
        let mut t = f + dt;
        for _ii in 1..=8i32 {
            let daux = CurveEval::derivative_at(the_curve, t);
            norm += daux.length();
            t += dt;
        }
        norm /= 9.0;
        dt = d1.length() / norm;
        if dt < 1.5 && dt > 0.75 {
            // The edge is within the average, keep it.
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
    cons.set(1, size as i32, the_point.x);
    cons.set(2, size as i32, the_point.y);
    cons.set(3, size as i32, the_point.z);
    if the_continuity >= 2 {
        cons.set(1, 3, d2.x);
        cons.set(2, 3, d2.y);
        cons.set(3, 3, d2.z);
    }
    if the_continuity >= 3 {
        cons.set(1, 4, d3.x);
        cons.set(2, 4, d3.y);
        cons.set(3, 4, d3.z);
    }
    // OCCT L1368.
    compute_lambda(&cons, &mat_coefs, l1, &mut lambda);

    // Construction in the Polynomial Basis (OCCT L1370-1394).
    let mut cont = vec![DVec3::ZERO; size as usize];
    cont[0] = p0;
    cont[1] = d1 * lambda;
    if the_continuity >= 2 {
        cont[2] = d2 * lambda.powi(2);
    }
    if the_continuity >= 3 {
        cont[3] = d3 * lambda.powi(3);
    }
    cont[size as usize - 1] = the_point;

    // NCollection_Array1<gp_Pnt> ExtrapPoles / ExtraCoeffs (OCCT L1384-1394).
    let mut extra_coeffs = vec![DVec3::ZERO; size as usize];
    for ii in 1..=size as usize {
        for jj in 1..=size as usize {
            extra_coeffs[jj - 1] += mat_coefs.get(ii as i32, jj as i32) * cont[ii - 1];
        }
    }

    // Conversion to the Bernstein Basis (OCCT L1396-1397):
    // PLib::CoefficientsPoles(ExtraCoeffs, NoWeights, ExtrapPoles, NoWeights).
    let extrap_poles = plib_coefficients_poles_dim3(&extra_coeffs);

    // OCCT L1399: Geom_BezierCurve Bezier = new Geom_BezierCurve(ExtrapPoles).
    let a_bezier = Curve3::Bezier(BezierCurve3 {
        control_points: extrap_poles.clone(),
        weights: vec![1.0; extrap_poles.len()],
    });

    // OCCT L1401-1404: Tol += ExtrapPoles(1).Distance(p0).
    let dist = extrap_poles[0].distance(p0);
    tol += dist;

    // Concatenation (OCCT L1406-1409).
    let ok = concat.add(&a_bezier, tol, the_after);
    if !ok {
        panic!("Standard_ConstructionError: ExtendCurveToPoint");
    }
    // Curve = Concat.BSplineCurve();
    *the_curve = Curve3::BSpline(
        concat
            .bspline_curve()
            .expect("GeomConvert_CompCurveToBSplineCurve::BSplineCurve (null)"),
    );
}

// ---------------------------------------------------------------------------
// OCCT GeomLib::To3d (GeomLib.cxx L559-675) — lift a 2D curve into the 3D
// frame of a gp_Ax2. The analytic overloads delegate to the kernel
// `ElCLib::To3d` set of `rcad_kernel::math::el`.
// ---------------------------------------------------------------------------

/// OCCT `GeomLib::To3d(const gp_Ax2& Position, const Handle(Geom2d_Curve)&)`
/// (GeomLib.cxx L559-675).
///
/// TrimmedCurve: recurse into the basis and re-trim (L565-573);
/// OffsetCurve: recurse and rebuild `Geom_OffsetCurve(CC, Offset,
/// Position.Direction())` (L574-581); Bezier / BSpline: map the poles through
/// `ElCLib::To3d` keeping the weights / knots / degree (L582-636); Line /
/// Circle / Ellipse / Parabola / Hyperbola: the `ElCLib::To3d` frame lift
/// (L637-671). Any other type raises `Standard_NotImplemented` in OCCT
/// (L672-674).
///
/// Architecture difference: the rcad `BSplineCurve2` carries no periodicity
/// flag, so the 3D `is_periodic` is derived from the knot structure (the same
/// predicate `rcad_kernel::math::bspl::bspline_is_periodic` uses for the OCCT
/// accessor).
pub fn to_3d(position: &rcad_kernel::math::gp::Ax2, curve2d: &rcad_kernel::geom::Curve2d) -> Option<Curve3> {
    use rcad_kernel::math::el::{
        elclib_to3d_circle, elclib_to3d_ellipse, elclib_to3d_hyperbola, elclib_to3d_line,
        elclib_to3d_parabola, elclib_to3d_pnt,
    };
    use rcad_kernel::geom::Curve2d as C2d;
    match curve2d {
        // L565-573: TrimmedCurve — recurse, then re-trim.
        C2d::Trimmed(tc) => {
            let cc = to_3d(position, &tc.curve)?;
            Some(Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3::new(
                cc,
                tc.t_min,
                tc.t_max,
            )))
        }
        // L574-581: OffsetCurve.
        C2d::Offset(co) => {
            let cc = to_3d(position, &co.basis)?;
            Some(Curve3::Offset(rcad_kernel::geom::OffsetCurve3 {
                basis: Box::new(cc),
                offset_distance: co.offset_distance,
                offset_dir: position.direction,
            }))
        }
        // L582-601: BezierCurve — pole mapping, weights kept.
        C2d::Bezier(b) => {
            let poles: Vec<DVec3> = b
                .control_points
                .iter()
                .map(|p| elclib_to3d_pnt(position, *p))
                .collect();
            Some(Curve3::Bezier(BezierCurve3 {
                control_points: poles,
                weights: b.weights.clone(),
            }))
        }
        // L603-636: BSplineCurve — pole mapping; degree / knots / weights are
        // carried by the rcad BSplineCurve2 storage (the flat knot vector
        // expands the OCCT multiplicities).
        C2d::BSpline(b) => {
            let poles: Vec<DVec3> = b
                .control_points
                .iter()
                .map(|p| elclib_to3d_pnt(position, *p))
                .collect();
            Some(Curve3::BSpline(rcad_kernel::geom::BSplineCurve3 {
                degree: b.degree,
                knots: b.knots.clone(),
                control_points: poles,
                weights: b.weights.clone(),
                is_periodic: rcad_kernel::math::bspl::bspline_is_periodic(&b.knots, b.degree),
            }))
        }
        // L637-643: Line2d.
        C2d::Line(l) => Some(Curve3::Line(elclib_to3d_line(
            position,
            l.origin,
            l.direction,
        ))),
        // L644-651: Circle2d.
        C2d::Circle(c) => Some(Curve3::Circle(elclib_to3d_circle(position, c))),
        // L652-659: Ellipse2d.
        C2d::Ellipse(e) => Some(Curve3::Ellipse(elclib_to3d_ellipse(position, e))),
        // L660-666: Parabola2d.
        C2d::Parabola(p) => Some(Curve3::Parabola(elclib_to3d_parabola(position, p))),
        // L667-671: Hyperbola2d.
        C2d::Hyperbola(h) => Some(Curve3::Hyperbola(elclib_to3d_hyperbola(position, h))),
        // L672-674: throw Standard_NotImplemented().
        _ => panic!("Standard_NotImplemented: GeomLib::To3d (GeomLib.cxx L672-674)"),
    }
}

// ---------------------------------------------------------------------------
// OCCT GeomLib_CurveOnSurfaceEvaluator (GeomLib.cxx L974-1050)
// ---------------------------------------------------------------------------

/// OCCT GeomLib_CurveOnSurfaceEvaluator (GeomLib.cxx L976-999) — the
/// `AdvApprox_EvaluatorFunction` over an `Adaptor3d_CurveOnSurface` restricted
/// to the current [First, Last].
///
/// The OCCT class holds `Adaptor3d_CurveOnSurface& CurveOnSurface` (a
/// reference; the rcad translation keeps the same borrow) plus the cached
/// `handle(Adaptor3d_Curve) TrimCurve` and the current `FirstParam` /
/// `LastParam`.
pub struct GeomLibCurveOnSurfaceEvaluator<'a> {
    /// OCCT: Adaptor3d_CurveOnSurface& CurveOnSurface.
    curve_on_surface: &'a CurveOnSurface,
    /// OCCT: double FirstParam.
    first_param: f64,
    /// OCCT: double LastParam.
    last_param: f64,
    /// OCCT: handle(Adaptor3d_Curve) TrimCurve.
    trim_curve: Option<Arc<dyn Adaptor3dCurve>>,
}

impl<'a> GeomLibCurveOnSurfaceEvaluator<'a> {
    /// OCCT ctor (GeomLib.cxx L979-986).
    pub fn new(curve_on_surface: &'a CurveOnSurface, the_first: f64, the_last: f64) -> Self {
        GeomLibCurveOnSurfaceEvaluator {
            curve_on_surface,
            first_param: the_first,
            last_param: the_last,
            trim_curve: None,
        }
    }
}

impl EvaluatorFunction for GeomLibCurveOnSurfaceEvaluator<'_> {
    /// OCCT GeomLib_CurveOnSurfaceEvaluator::Evaluate (GeomLib.cxx L1001-1050).
    fn evaluate(
        &mut self,
        start_end: &[f64; 2],
        parameter: f64,
        derivative_request: i32,
        result: &mut [f64],
    ) -> i32 {
        // Handle left / right positioning (cxx L1007-1014).
        if start_end[0] != self.first_param || start_end[1] != self.last_param {
            // OCCT: TrimCurve = CurveOnSurface.Trim(DebutFin[0], DebutFin[1],
            // Precision::PConfusion()).
            self.trim_curve = Some(self.curve_on_surface.trim(
                start_end[0],
                start_end[1],
                p_confusion(),
            ));
            self.first_param = start_end[0];
            self.last_param = start_end[1];
        }

        // Positioning (cxx L1016-1043).
        let trim_curve = self.trim_curve.as_ref().expect(
            "Standard_NullObject: GeomLib_CurveOnSurfaceEvaluator TrimCurve",
        );
        if derivative_request == 0 {
            let point = trim_curve.value(parameter);
            for ii in 0..3 {
                result[ii] = point[ii];
            }
        }
        if derivative_request == 1 {
            let (_point, vector) = trim_curve.d1(parameter);
            for ii in 0..3 {
                result[ii] = vector[ii];
            }
        }
        if derivative_request == 2 {
            // OCCT: TrimCurve->D2((*Parameter), Point, VecBis, Vector) — the
            // second derivative argument order of Adaptor3d_Curve::D2.
            let (_point, _vec_bis, vector) = trim_curve.d2(parameter);
            for ii in 0..3 {
                result[ii] = vector[ii];
            }
        }
        // ReturnCode[0] = 0 (cxx L1045).
        0
    }
}

/// OCCT GeomLib::BuildCurve3d(Tolerance, Curve, FirstParameter, LastParameter,
/// NewCurvePtr, MaxDeviation, AverageDeviation, Continuity, MaxDegree,
/// MaxSegment) (GeomLib.cxx L1051-1163).
///
/// `Curve` is rcad's kernel `CurveOnSurface` — the real
/// `Adaptor3d_CurveOnSurface`.  `new_curve_ptr` is the OCCT `NewCurvePtr`
/// out-parameter: a null handle on every OCCT failure path (the failed
/// down-casts, the failed isoline construction, and the no-result
/// approximation).
#[allow(clippy::too_many_arguments)]
pub fn build_curve3d(
    tolerance: f64,
    curve: &CurveOnSurface,
    first_parameter: f64,
    last_parameter: f64,
    new_curve_ptr: &mut Option<Curve3>,
    max_deviation: &mut f64,
    average_deviation: &mut f64,
    continuity: GeomAbsShape,
    max_degree: i32,
    max_segment: i32,
) {
    // OCCT L1064-1065.
    *max_deviation = 0.0e0;
    *average_deviation = 0.0e0;

    // OCCT L1066-1070: the two down-casts.  The rcad bridges model the OCCT
    // null-handle outcome of a failed down-cast as `None`.
    let geom_adaptor_surface_ptr: Option<&Surface3> = curve.get_surface().kernel_surface();
    let geom_adaptor_curve_ptr: Option<&Curve2d> = curve.get_curve().kernel_curve2d();

    if let (Some(geom2d_curve), Some(geom_surface)) =
        (geom_adaptor_curve_ptr, geom_adaptor_surface_ptr)
    {
        // OCCT L1076-1090: the Geom_RectangularTrimmedSurface unwrap and the
        // Geom_Plane detection.
        let p: Option<&rcad_kernel::geom::Plane> = match geom_surface {
            Surface3::Trimmed(rt) => match rt.basis.as_ref() {
                Surface3::Plane(a_plane) => Some(a_plane),
                _ => None,
            },
            Surface3::Plane(a_plane) => Some(a_plane),
            _ => None,
        };

        if let Some(a_plane) = p {
            // OCCT L1094-1098: compute the 3d curve.
            let axes = rcad_kernel::math::gp::Ax2::new(
                a_plane.origin,
                a_plane.normal,
                a_plane.u_dir,
            );
            *new_curve_ptr = to_3d(&axes, geom2d_curve);
            return;
        }

        // OCCT L1100-1101: TrimmedC2D = geom_adaptor_curve_ptr->Trim(
        // FirstParameter, LastParameter, Precision::PConfusion()).
        let trimmed_c2d: Curve2dHandle =
            curve
                .get_curve()
                .trim(first_parameter, last_parameter, p_confusion());

        let mut is_u = false;
        let mut a_param = 0.0f64;
        let mut is_forward = false;
        if is_iso_line(trimmed_c2d.as_ref(), &mut is_u, &mut a_param, &mut is_forward) {
            // OCCT L1105-1113: NewCurvePtr = buildC3dOnIsoLine(TrimmedC2D,
            // geom_adaptor_surface_ptr, FirstParameter, LastParameter,
            // Tolerance, isU, aParam, isForward).
            let built = build_c3d_on_iso_line(
                trimmed_c2d.as_ref(),
                curve.get_surface().as_ref(),
                first_parameter,
                last_parameter,
                tolerance,
                is_u,
                a_param,
                is_forward,
            );
            if let Some(a_c3d) = built {
                *new_curve_ptr = Some(a_c3d);
                return;
            }
        }
    }

    //
    // Entree
    //
    // OCCT L1118-1121: Tolerance1DPtr / Tolerance2DPtr stay null, Tolerance3DPtr
    // is the single 3D tolerance.
    let tolerance3d: Vec<f64> = vec![tolerance];

    // Search for discontinuities (OCCT L1123-1130).
    let _nb_interval_c2 = curve.nb_intervals(GeomAbsShape::C2);
    let param_de_decoupe_c2: Vec<f64> = curve.intervals(GeomAbsShape::C2);

    let _nb_interval_c3 = curve.nb_intervals(GeomAbsShape::C3);
    let param_de_decoupe_c3: Vec<f64> = curve.intervals(GeomAbsShape::C3);

    // Note extension of the parametric range.
    // To force Trim on first evaluator call.
    // OCCT L1132-1134.
    let mut ev = GeomLibCurveOnSurfaceEvaluator::new(
        curve,
        first_parameter - 1.0,
        last_parameter + 1.0,
    );

    // Approximation with preferential cutting (OCCT L1136-1149).
    let preferentiel =
        PrefAndRec::with_default_weight(&param_de_decoupe_c2, &param_de_decoupe_c3);
    let an_approximator = ApproxAFunction::with_cut_tool(
        0,
        0,
        1,
        None,
        None,
        Some(&tolerance3d),
        first_parameter,
        last_parameter,
        continuity,
        max_degree,
        max_segment,
        &mut ev,
        &preferentiel,
    );

    if an_approximator.has_result() {
        // OCCT L1153-1161.
        let a_curve_builder = GeomLibMakeCurvefromApprox::new(&an_approximator);

        let a_curve_ptr = a_curve_builder.curve(1);
        // Return the approximation results.
        *max_deviation = an_approximator.max_error_at(3, 1);
        *average_deviation = an_approximator.average_error_at(3, 1);
        *new_curve_ptr = a_curve_ptr.map(Curve3::BSpline);
    }
}
