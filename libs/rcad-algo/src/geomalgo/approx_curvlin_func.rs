//! OCCT Approx_CurvlinFunc (ModelingData/TKGeomBase/Approx —
//! Approx_CurvlinFunc.hxx L27-141, Approx_CurvlinFunc.cxx L44-793).
//!
//! 1:1 translation of the curvilinear-parametrization function class consumed
//! by Approx_CurvilinearParameter (the BRepLib::SameParameter IsBad branch).
//!
//! rcad encoding of the OCCT handles (architecture differences, annotated):
//! - `handle(Adaptor3d_Curve) myC3D` -> [`CurveHandle`]
//!   (`Arc<dyn Adaptor3dCurve>`),
//! - `handle(Adaptor2d_Curve2d) myC2D1/myC2D2` -> [`Curve2dHandle`],
//! - `handle(Adaptor3d_Surface) mySurf1/mySurf2` -> [`SurfaceHandle`],
//! - `handle(NCollection_HArray1<double>) myUi_1/mySi_1/myUi_2/mySi_2` ->
//!   `Option<Vec<f64>>` (the OCCT arrays are 0-based `(0, NbIntC3 * NbInt)`;
//!   the Vec index is the same 0-based ordinal),
//! - the OCCT `const_cast` writes of myPrevS/myPrevU inside the const
//!   `Init`/`GetUParameter` map to `&mut self` receivers (the cache fields
//!   are the only mutations).
//!
//! The file additionally re-hosts the GCPnts/CPnts machinery the OCCT
//! GetUParameter relies on and that the
//! [`gcpnts_abscissa_point`](super::gcpnts_abscissa_point) module does not
//! carry yet (its header records the instance machinery as not delivered):
//! - `math_GaussSingleIntegration` (math_GaussSingleIntegration.cxx L64-155),
//! - `CPnts_MyGaussFunction` / `CPnts_MyRootFunction` (CPnts_MyRootFunction.hxx
//!   L38-62, CPnts_MyRootFunction.cxx L30-96),
//! - `CPnts_AbscissaPoint` instance machinery (CPnts_AbscissaPoint.cxx
//!   L270-474: Init / Perform / AdvPerform),
//! - `GCPnts_AbscissaPoint` 5-argument constructor through `AdvCompute`
//!   (GCPnts_AbscissaPoint.cxx L26-65 computeType, L166-297 AdvCompute,
//!   L533-553 the (C, Abscissa, U0, Ui, Tol) constructor).

use std::sync::Arc;

use glam::DVec3;

use rcad_kernel::base::gcpnts::gcpnts_curve::GCPntsCurve;
use rcad_kernel::base::proj_lib::adaptor::{
    Adaptor2dCurve2d, Adaptor3dCurve, Curve2dHandle, CurveHandle, SurfaceHandle,
};
use rcad_kernel::base::proj_lib::CurveType;
use rcad_kernel::base::geom_lib::fuse_intervals;
use rcad_kernel::math::gauss_points::{gauss_points, gauss_points_max, gauss_weights};
use rcad_kernel::math::root::{FunctionValue, FunctionWithDerivative};
use rcad_kernel::math::VecD;
use rcad_kernel::math::{bspl_lib, GeomAbsShape};

use crate::geomalgo::gcpnts_abscissa_point::gcpnts_length_3d_range_tol;

// =========================================================================
// GCPnts/CPnts re-hosts (see the module header for the hosting rationale).
// =========================================================================

/// OCCT math_GaussSingleIntegration::Perform
/// (math_GaussSingleIntegration.cxx L100-155) — the Gauss quadrature of
/// `value` over [lower, upper] with `order` points.
fn gauss_single_integration_perform(
    value: &mut dyn FnMut(f64) -> f64,
    lower: f64,
    upper: f64,
    order: usize,
) -> Option<f64> {
    let the_order = order.min(gauss_points_max());
    let mut gauss_p = VecD::new(the_order + 1);
    let mut gauss_w = VecD::new(the_order + 1);
    // math::GaussPoints(Order, GaussP); math::GaussWeights(Order, GaussW).
    gauss_points(the_order, &mut gauss_p);
    gauss_weights(the_order, &mut gauss_w);

    // Changement de variable pour la mise a l'echelle [Lower, Upper].
    let xm = 0.5 * (upper + lower);
    let xr = 0.5 * (upper - lower);
    let mut val = 0.0f64;

    let ind = the_order / 2;
    let ind1 = (the_order + 1) / 2;
    if ind1 > ind {
        // odder case
        val = value(xm);
        val *= gauss_w.get(ind1);
    }
    // Sommation sur tous les points de Gauss: avec utilisation de la symetrie.
    for j in 1..=ind {
        let dx = xr * gauss_p.get(j);
        let f1 = value(xm - dx);
        let f2 = value(xm + dx);
        // Multiplication par les poids de Gauss.
        let ft = f1 + f2;
        val += gauss_w.get(j) * ft;
    }
    // Mise a l'echelle de l'intervalle [Lower, Upper]
    val *= xr;
    Some(val)
}

/// OCCT math_GaussSingleIntegration(F, Lower, Upper, Order, Tol)
/// (math_GaussSingleIntegration.cxx L64-98) — the adaptive subdivision loop.
fn gauss_single_integration_tol(
    value: &mut dyn FnMut(f64) -> f64,
    lower: f64,
    upper: f64,
    order: usize,
    tol: f64,
) -> Option<f64> {
    let the_order = order.min(gauss_points_max());
    let iter_max = 13; // Max number of iteration
    let mut n_iter = 1; // current number of iteration
    let mut nb_interval = 1; // current number of subintervals

    let mut len = gauss_single_integration_perform(value, lower, upper, the_order)?;
    loop {
        let old_len = len;
        len = 0.0;
        nb_interval *= 2;
        let du = (upper - lower) / nb_interval as f64;
        for i in 1..=nb_interval {
            let v = gauss_single_integration_perform(
                value,
                lower + (i - 1) as f64 * du,
                lower + i as f64 * du,
                the_order,
            )?;
            len += v;
        }
        n_iter += 1;
        if !((old_len - len).abs() > tol && n_iter <= iter_max) {
            break;
        }
    }
    Some(len)
}

/// The OCCT CPnts_RealFunction callback — the magnitude integrand over an
/// `Adaptor3d_Curve` (CPnts_AbscissaPoint.cxx L41-47 f3d).
type CPntsRealFunction<'a> = Box<dyn Fn(f64) -> f64 + 'a>;

/// OCCT CPnts_MyGaussFunction (CPnts_MyGaussFunction.hxx L32-46) — the
/// math_Function adapter evaluating the integrand at X.
struct CPntsMyGaussFunction<'a> {
    /// OCCT: CPnts_RealFunction myFunction + void* myData — the rcad
    /// closure carries both (the data pointer is the captured curve).
    my_function: Option<CPntsRealFunction<'a>>,
}

impl<'a> CPntsMyGaussFunction<'a> {
    /// OCCT CPnts_MyGaussFunction() (CPnts_MyGaussFunction.lxx L17-20).
    fn new() -> Self {
        CPntsMyGaussFunction { my_function: None }
    }

    /// OCCT CPnts_MyGaussFunction::Init(F, D) (CPnts_MyGaussFunction.cxx).
    fn init(&mut self, the_f: CPntsRealFunction<'a>) {
        self.my_function = Some(the_f);
    }

    /// OCCT CPnts_MyGaussFunction::Value(X, F) (CPnts_MyGaussFunction.cxx).
    fn value(&self, the_x: f64) -> Option<f64> {
        match &self.my_function {
            Some(the_f) => Some(the_f(the_x)),
            None => None,
        }
    }
}

/// OCCT CPnts_MyRootFunction (CPnts_MyRootFunction.hxx L38-62) — the
/// math_FunctionWithDerivative solving Integral(X0,X,F(X)) = L.
struct CPntsMyRootFunction<'a> {
    /// OCCT: CPnts_MyGaussFunction myFunction.
    my_function: CPntsMyGaussFunction<'a>,
    /// OCCT: double myX0.
    my_x0: f64,
    /// OCCT: double myL.
    my_l: f64,
    /// OCCT: int myOrder.
    my_order: usize,
    /// OCCT: double myTol.
    my_tol: f64,
}

impl<'a> CPntsMyRootFunction<'a> {
    /// OCCT CPnts_MyRootFunction() (CPnts_MyRootFunction.lxx L17-22).
    fn new() -> Self {
        CPntsMyRootFunction {
            my_function: CPntsMyGaussFunction::new(),
            my_x0: 0.0,
            my_l: 0.0,
            my_order: 0,
            my_tol: 0.0,
        }
    }

    /// OCCT CPnts_MyRootFunction::Init(F, D, Order)
    /// (CPnts_MyRootFunction.cxx L30-34).
    fn init_function(&mut self, the_f: CPntsRealFunction<'a>, the_order: usize) {
        self.my_function.init(the_f);
        self.my_order = the_order;
    }

    /// OCCT CPnts_MyRootFunction::Init(X0, L) (CPnts_MyRootFunction.cxx
    /// L36-41) — myTol = -1 suppresses the tolerance.
    fn init_equation(&mut self, the_x0: f64, the_l: f64) {
        self.my_x0 = the_x0;
        self.my_l = the_l;
        self.my_tol = -1.0;
    }

    /// OCCT CPnts_MyRootFunction::Init(X0, L, Tol) (CPnts_MyRootFunction.cxx
    /// L43-48).
    #[allow(dead_code)]
    fn init_equation_tol(&mut self, the_x0: f64, the_l: f64, the_tol: f64) {
        self.my_x0 = the_x0;
        self.my_l = the_l;
        self.my_tol = the_tol;
    }
}

impl FunctionValue for CPntsMyRootFunction<'_> {
    /// OCCT CPnts_MyRootFunction::Value(X, F) (CPnts_MyRootFunction.cxx
    /// L50-72) — F = Integral(X0, X) - L.
    fn value(&mut self, the_x: f64) -> Option<f64> {
        let the_length = if self.my_tol <= 0.0 {
            gauss_single_integration_perform(
                &mut |t| self.my_function.value(t).unwrap_or(f64::NAN),
                self.my_x0,
                the_x,
                self.my_order,
            )
        } else {
            gauss_single_integration_tol(
                &mut |t| self.my_function.value(t).unwrap_or(f64::NAN),
                self.my_x0,
                the_x,
                self.my_order,
                self.my_tol,
            )
        };
        the_length.map(|v| v - self.my_l)
    }
}

impl FunctionWithDerivative for CPntsMyRootFunction<'_> {
    /// (the math_FunctionWithDerivative extension: Derivative + Values; the
    /// Value arm lives in the FunctionValue impl above.)
    /// OCCT CPnts_MyRootFunction::Derivative(X, Df) (CPnts_MyRootFunction.cxx
    /// L74-78) — Df = F(X).
    fn derivative(&mut self, the_x: f64) -> Option<f64> {
        self.my_function.value(the_x)
    }

    /// OCCT CPnts_MyRootFunction::Values(X, F, Df) (CPnts_MyRootFunction.cxx
    /// L80-105).
    fn values(&mut self, the_x: f64) -> Option<(f64, f64)> {
        let the_f = FunctionValue::value(self, the_x)?;
        let the_df = self.my_function.value(the_x)?;
        Some((the_f, the_df))
    }
}

/// OCCT CPnts_AbscissaPoint::order(C) — the Adaptor3d_Curve flavor
/// (CPnts_AbscissaPoint.cxx L57-77); no clamp in OCCT (the degree-0 /
/// single-pole degenerate cases cannot reach here through the callers).
fn cpnts_order_3d(the_c: &dyn GCPntsCurve) -> usize {
    match the_c.get_type() {
        CurveType::Line => 2,
        CurveType::Parabola => 5,
        CurveType::Bezier => (24).min(2 * the_c.curve_degree()) as usize,
        CurveType::BSpline => (24).min(2 * the_c.nb_poles() - 1) as usize,
        _ => 10,
    }
}

/// OCCT GCPnts_AbscissaType.hxx — compute/length classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GCPntsAbscissaType {
    AbsComposite,
    LengthParametrized,
    Parametrized,
}

/// OCCT static computeType (GCPnts_AbscissaPoint.cxx L26-65) — the type and
/// the length ratio if GCPnts_LengthParametrized.
fn gcpnts_compute_type(the_c: &dyn GCPntsCurve, the_ratio: &mut f64) -> GCPntsAbscissaType {
    if the_c.nb_intervals_cn() > 1 {
        return GCPntsAbscissaType::AbsComposite;
    }

    match the_c.get_type() {
        CurveType::Line => {
            *the_ratio = 1.0;
            GCPntsAbscissaType::LengthParametrized
        }
        CurveType::Circle => {
            *the_ratio = the_c.circle_radius();
            GCPntsAbscissaType::LengthParametrized
        }
        CurveType::Bezier => {
            if the_c.nb_poles() == 2 && !the_c.is_rational() {
                *the_ratio = the_c.dn1(0.0).length();
                GCPntsAbscissaType::LengthParametrized
            } else {
                GCPntsAbscissaType::Parametrized
            }
        }
        CurveType::BSpline => {
            if the_c.nb_poles() == 2 && !the_c.is_rational() {
                *the_ratio = the_c.dn1(the_c.first_parameter()).length();
                GCPntsAbscissaType::LengthParametrized
            } else {
                GCPntsAbscissaType::Parametrized
            }
        }
        _ => GCPntsAbscissaType::Parametrized,
    }
}

/// OCCT CPnts_AbscissaPoint (CPnts_AbscissaPoint.cxx L211-474) — the instance
/// machinery (myDone / myL / myParam / myUMin / myUMax + MyRootFunction).
struct CPntsAbscissaPoint<'a> {
    /// OCCT: CPnts_MyRootFunction myF.
    my_f: CPntsMyRootFunction<'a>,
    /// OCCT: bool myDone.
    my_done: bool,
    /// OCCT: double myL.
    my_l: f64,
    /// OCCT: double myParam.
    my_param: f64,
    /// OCCT: double myUMin.
    my_umin: f64,
    /// OCCT: double myUMax.
    my_umax: f64,
}

impl<'a> CPntsAbscissaPoint<'a> {
    /// OCCT CPnts_AbscissaPoint() (CPnts_AbscissaPoint.cxx L211-218).
    fn new() -> Self {
        CPntsAbscissaPoint {
            my_f: CPntsMyRootFunction::new(),
            my_done: false,
            my_l: 0.0,
            my_param: 0.0,
            my_umin: 0.0,
            my_umax: 0.0,
        }
    }

    /// OCCT CPnts_AbscissaPoint::SetParameter(U) — the accessor used by the
    /// Compute drivers (the myParam write, CPnts_AbscissaPoint.hxx inline).
    fn set_parameter(&mut self, the_u: f64) {
        self.my_param = the_u;
        self.my_done = true;
    }

    /// OCCT CPnts_AbscissaPoint::Init(C, U1, U2, Tol) — the Adaptor3d_Curve
    /// flavor (CPnts_AbscissaPoint.cxx L336-351).
    fn init(
        &mut self,
        the_c: &'a dyn Adaptor3dCurve,
        the_u1: f64,
        the_u2: f64,
        the_tol: f64,
    ) {
        // OCCT L342-343: myF.Init(f3d, (void*)&C, order(C)).
        let a_order = cpnts_order_3d(&AdaptorAsGCPnts(the_c));
        self.my_f.init_function(Box::new(move |x| the_c.d1(x).1.length()), a_order);
        // OCCT L344: myL = CPnts_AbscissaPoint::Length(C, U1, U2, Tol).
        self.my_l = gcpnts_length_3d_range_tol(
            &AdaptorAsGCPnts(the_c),
            the_u1,
            the_u2,
            the_tol,
        );
        // OCCT L345-350.
        self.my_umin = the_u1.min(the_u2);
        self.my_umax = the_u1.max(the_u2);
        let a_du = self.my_umax - self.my_umin;
        self.my_umin -= a_du;
        self.my_umax += a_du;
    }

    /// OCCT CPnts_AbscissaPoint::AdvPerform(Abscissa, U0, Ui, Resolution)
    /// (CPnts_AbscissaPoint.cxx L436-474) — rbv's tolerance management.
    fn adv_perform(&mut self, the_abscissa: f64, the_u0: f64, the_ui: f64, the_resolution: f64) {
        if self.my_l < rcad_kernel::precision::CONFUSION {
            // leave less violently :
            self.my_done = true;
            self.my_param = the_u0;
        } else {
            self.my_done = false;
            // OCCT L453: myF.Init(U0, Abscissa, Resolution / 10).
            self.my_f.init_equation_tol(the_u0, the_abscissa, the_resolution / 10.0);
            // OCCT L455: math_FunctionRoot Solution(myF, Ui, Resolution,
            //               myUMin, myUMax) — the bounded form with the
            // math_FunctionRoot.hxx L46 default NbIterations = 100.
            let a_solution = rcad_kernel::math::newton_function_root::NewtonFunctionRoot::new_bounded(
                &mut self.my_f,
                the_ui,
                the_resolution,
                self.my_umin,
                self.my_umax,
                100,
            );
            // OCCT L468-472.
            if a_solution.is_done() {
                self.my_done = true;
                self.my_param = a_solution.root();
            }
        }
    }

    /// OCCT CPnts_AbscissaPoint::Parameter().
    fn parameter(&self) -> f64 {
        self.my_param
    }
}

/// Architecture glue (no direct OCCT counterpart; the same encoding as the
/// kernel `CurveToolAsGCPnts` bridge): the `Adaptor3d_Curve` flavor of the
/// GCPnts template interface served over a `&dyn Adaptor3dCurve` (OCCT
/// instantiates the GCPnts templates with `TheCurve = Adaptor3d_Curve`).
struct AdaptorAsGCPnts<'a>(&'a dyn Adaptor3dCurve);

impl GCPntsCurve for AdaptorAsGCPnts<'_> {
    fn first_parameter(&self) -> f64 {
        self.0.first_parameter()
    }

    fn last_parameter(&self) -> f64 {
        self.0.last_parameter()
    }

    fn d0(&self, u: f64) -> DVec3 {
        self.0.value(u)
    }

    fn d1(&self, u: f64) -> (DVec3, DVec3) {
        self.0.d1(u)
    }

    fn d2(&self, u: f64) -> (DVec3, DVec3, DVec3) {
        self.0.d2(u)
    }

    /// OCCT GetType() — the erased handle answers through the trait
    /// `curve_type` accessor (the GeomCurveAdaptor override reads the
    /// loaded myTypeCurve; the trait default is the OCCT base
    /// GeomAbs_OtherCurve).
    fn get_type(&self) -> CurveType {
        self.0.curve_type()
    }

    fn nb_intervals_cn(&self) -> usize {
        self.0.nb_intervals(GeomAbsShape::CN)
    }

    fn intervals_cn(&self) -> Vec<f64> {
        self.0.intervals(GeomAbsShape::CN)
    }

    fn circle_radius(&self) -> f64 {
        panic!("GCPntsCurve: Circle() on a non-circle curve");
    }

    fn curve_degree(&self) -> i32 {
        self.0.degree() as i32
    }

    fn nb_poles(&self) -> i32 {
        self.0.nb_poles() as i32
    }

    fn is_rational(&self) -> bool {
        self.0.is_rational()
    }

    fn dn1(&self, u: f64) -> DVec3 {
        self.0.dn(u, 1)
    }
}

/// The GCPnts_AbscissaPoint(C, Abscissa, U0, Ui, Tol) 5-argument constructor
/// (GCPnts_AbscissaPoint.cxx L533-541) — AdvCompute + Parameter() as one
/// step (the rcad value encoding of the OCCT constructor object).
fn gcpnts_abscissa_point_adv<'a>(
    the_c: &'a dyn Adaptor3dCurve,
    the_abscissa: f64,
    the_u0: f64,
    the_ui: f64,
    the_tol: f64,
) -> f64 {
    let mut a_computer = CPntsAbscissaPoint::new();
    let mut an_abscis = the_abscissa;
    let mut a_uu0 = the_u0;
    let mut a_uui = the_ui;
    gcpnts_adv_compute(the_tol, the_c, &mut an_abscis, &mut a_uu0, &mut a_uui, &mut a_computer);
    a_computer.parameter()
}

/// OCCT static AdvCompute (GCPnts_AbscissaPoint.cxx L166-297) — performs
/// more appropriate tolerance management (introduced by rbv for curvilinear
/// parametrization).
fn gcpnts_adv_compute<'a>(
    the_epsilon: f64,
    the_c: &'a dyn Adaptor3dCurve,
    the_abscis: &mut f64,
    the_u0: &mut f64,
    the_ui: &mut f64,
    the_computer: &mut CPntsAbscissaPoint<'a>,
) {
    let a_curve = AdaptorAsGCPnts(the_c);
    let mut a_ratio = 1.0f64;
    let a_type = gcpnts_compute_type(&a_curve, &mut a_ratio);
    match a_type {
        GCPntsAbscissaType::LengthParametrized => {
            the_computer.set_parameter(*the_u0 + *the_abscis / a_ratio);
        }
        GCPntsAbscissaType::Parametrized => {
            // theComputer.Init(theC, theEPSILON) — rbv's modification
            // (GCPnts_AbscissaPoint.cxx L183).
            the_computer.init(the_c, the_c.first_parameter(), the_c.last_parameter(), the_epsilon);
            the_computer.adv_perform(*the_abscis, *the_u0, *the_ui, the_epsilon);
        }
        GCPntsAbscissaType::AbsComposite => {
            let a_nb_intervals = the_c.nb_intervals(GeomAbsShape::CN);
            let a_ti = the_c.intervals(GeomAbsShape::CN);
            let mut a_l = 0.0f64;
            let mut a_sign = 1.0f64;
            let mut an_index = 1i32;
            bspl_lib::hunt(&a_ti, *the_u0, &mut an_index);

            let mut a_direction = 1i32;
            if *the_abscis < 0.0 {
                a_direction = 0;
                *the_abscis = -*the_abscis;
                a_sign = -1.0;
            }

            if an_index == 0 && a_direction > 0 {
                a_l = gcpnts_length_3d_range_tol(
                    &a_curve,
                    *the_u0,
                    at(&a_ti, an_index + a_direction),
                    the_epsilon,
                );
                if (a_l - *the_abscis).abs() <= the_epsilon {
                    the_computer.set_parameter(at(&a_ti, an_index + a_direction));
                    return;
                }

                if a_l > *the_abscis {
                    if *the_ui > at(&a_ti, an_index + 1) {
                        *the_ui = (*the_abscis / a_l) * (at(&a_ti, an_index + 1) - *the_u0);
                        *the_ui = *the_u0 + *the_ui;
                    }
                    the_computer.init(the_c, *the_u0, at(&a_ti, an_index + 1), the_epsilon);
                    the_computer.adv_perform(a_sign * *the_abscis, *the_u0, *the_ui, the_epsilon);
                    return;
                } else {
                    *the_u0 = at(&a_ti, an_index + a_direction);
                    *the_abscis -= a_l;
                }
                an_index += 1;
            }

            while an_index >= 1 && an_index <= a_nb_intervals as i32 {
                a_l = gcpnts_length_3d_range_tol(
                    &a_curve,
                    *the_u0,
                    at(&a_ti, an_index + a_direction),
                    the_epsilon,
                );
                if (a_l - *the_abscis).abs() <= rcad_kernel::precision::PCONFUSION {
                    the_computer.set_parameter(at(&a_ti, an_index + a_direction));
                    return;
                }

                if a_l > *the_abscis {
                    if *the_ui < at(&a_ti, an_index) || *the_ui > at(&a_ti, an_index + 1) {
                        *the_ui = (*the_abscis / a_l) * (at(&a_ti, an_index + 1) - *the_u0);
                        if a_direction != 0 {
                            *the_ui = *the_u0 + *the_ui;
                        } else {
                            *the_ui = *the_u0 - *the_ui;
                        }
                    }
                    the_computer.init(
                        the_c,
                        at(&a_ti, an_index),
                        at(&a_ti, an_index + 1),
                        the_epsilon,
                    );
                    the_computer.adv_perform(a_sign * *the_abscis, *the_u0, *the_ui, the_epsilon);
                    return;
                } else {
                    *the_u0 = at(&a_ti, an_index + a_direction);
                    *the_abscis -= a_l;
                }
                if a_direction != 0 {
                    an_index += 1;
                } else {
                    an_index -= 1;
                }
            }

            // Push a little bit outside the limits (hairy !!!)
            let is_non_periodic = !the_c.is_periodic();
            *the_ui = *the_u0 + a_sign * 0.1;
            let mut a_u1 = *the_u0 + a_sign * 0.2;
            if is_non_periodic {
                if a_sign > 0.0 {
                    *the_ui = the_ui.min(the_c.last_parameter());
                    a_u1 = a_u1.min(the_c.last_parameter());
                } else {
                    *the_ui = the_ui.max(the_c.first_parameter());
                    a_u1 = a_u1.max(the_c.first_parameter());
                }
            }

            the_computer.init(the_c, *the_u0, a_u1, the_epsilon);
            the_computer.adv_perform(a_sign * *the_abscis, *the_u0, *the_ui, the_epsilon);
        }
    }
}

/// OCCT 1-based array read helper (`NCollection_Array1::Value(i)`).
fn at(v: &[f64], i: i32) -> f64 {
    v[(i - 1) as usize]
}

// =========================================================================
// Approx_CurvlinFunc file statics (Approx_CurvlinFunc.cxx L44-108).
// =========================================================================

/// OCCT static cubic (Approx_CurvlinFunc.cxx L44-60) — the Newton divided
/// differences of order 3 over (Xi, Yi) evaluated at X.
fn cubic(the_x: f64, the_xi: &[f64; 4], the_yi: &[f64; 4]) -> f64 {
    let i1 = (the_yi[0] - the_yi[1]) / (the_xi[0] - the_xi[1]);
    let i2 = (the_yi[1] - the_yi[2]) / (the_xi[1] - the_xi[2]);
    let i3 = (the_yi[2] - the_yi[3]) / (the_xi[2] - the_xi[3]);

    let i21 = (i1 - i2) / (the_xi[0] - the_xi[2]);
    let i22 = (i2 - i3) / (the_xi[1] - the_xi[3]);

    let i31 = (i21 - i22) / (the_xi[0] - the_xi[3]);

    the_yi[0] + (the_x - the_xi[0]) * (i1 + (the_x - the_xi[1]) * (i21 + (the_x - the_xi[2]) * i31))
}

/// OCCT static findfourpoints (Approx_CurvlinFunc.cxx L62-108) — picks the
/// four (S, U) knots around NInterval and tries to insert (prevS, prevU).
fn findfourpoints(
    _s: f64,
    mut n_interval: i32,
    the_si: &[f64],
    the_ui: &[f64],
    the_prev_s: f64,
    the_prev_u: f64,
    the_xi: &mut [f64; 4],
    the_yi: &mut [f64; 4],
) {
    let nb_int = the_si.len() as i32 - 1;
    if nb_int < 3 {
        // throw Standard_ConstructionError("Approx_CurvlinFunc::GetUParameter").
        panic!("Standard_ConstructionError: Approx_CurvlinFunc::GetUParameter");
    }

    if n_interval < 1 {
        n_interval = 1;
    } else if n_interval > nb_int - 2 {
        n_interval = nb_int - 2;
    }

    for a_i in 0..4i32 {
        the_xi[a_i as usize] = the_si[(n_interval - 1 + a_i) as usize];
        the_yi[a_i as usize] = the_ui[(n_interval - 1 + a_i) as usize];
    }
    // try to insert (S, U)
    for i in 0..3 {
        if the_xi[i] < the_prev_s && the_prev_s < the_xi[i + 1] {
            for j in 0..i {
                the_xi[j] = the_xi[j + 1];
                the_yi[j] = the_yi[j + 1];
            }
            the_xi[i] = the_prev_s;
            the_yi[i] = the_prev_u;
            break;
        }
    }
}

// =========================================================================
// Approx_CurvlinFunc (hxx L28-141).
// =========================================================================

/// OCCT Approx_CurvlinFunc (hxx L28-141).
pub struct ApproxCurvlinFunc {
    /// hxx L117: handle(Adaptor3d_Curve) myC3D.
    my_c3d: Option<CurveHandle>,
    /// hxx L118: handle(Adaptor2d_Curve2d) myC2D1.
    my_c2d1: Option<Curve2dHandle>,
    /// hxx L119: handle(Adaptor2d_Curve2d) myC2D2.
    my_c2d2: Option<Curve2dHandle>,
    /// hxx L120: handle(Adaptor3d_Surface) mySurf1.
    my_surf1: Option<SurfaceHandle>,
    /// hxx L121: handle(Adaptor3d_Surface) mySurf2.
    my_surf2: Option<SurfaceHandle>,
    /// hxx L122: int myCase.
    my_case: i32,
    /// hxx L123: double myFirstS.
    my_first_s: f64,
    /// hxx L124: double myLastS.
    my_last_s: f64,
    /// hxx L125: double myFirstU1.
    my_first_u1: f64,
    /// hxx L126: double myLastU1.
    my_last_u1: f64,
    /// hxx L127: double myFirstU2.
    my_first_u2: f64,
    /// hxx L128: double myLastU2.
    my_last_u2: f64,
    /// hxx L129: double myLength.
    my_length: f64,
    /// hxx L130: double myLength1.
    my_length1: f64,
    /// hxx L131: double myLength2.
    my_length2: f64,
    /// hxx L132: double myTolLen.
    my_tol_len: f64,
    /// hxx L133: double myPrevS.
    my_prev_s: f64,
    /// hxx L134: double myPrevU.
    my_prev_u: f64,
    /// hxx L135: handle(NCollection_HArray1<double>) myUi_1.
    my_ui_1: Option<Vec<f64>>,
    /// hxx L136: handle(NCollection_HArray1<double>) mySi_1.
    my_si_1: Option<Vec<f64>>,
    /// hxx L137: handle(NCollection_HArray1<double>) myUi_2.
    my_ui_2: Option<Vec<f64>>,
    /// hxx L138: handle(NCollection_HArray1<double>) mySi_2.
    my_si_2: Option<Vec<f64>>,
}

impl ApproxCurvlinFunc {
    /// OCCT Approx_CurvlinFunc(C, Tol) (cxx L127-137) — the case-1
    /// constructor.
    pub fn new_curve(the_c: CurveHandle, the_tol: f64) -> Self {
        let mut this = ApproxCurvlinFunc {
            my_c3d: Some(the_c),
            my_c2d1: None,
            my_c2d2: None,
            my_surf1: None,
            my_surf2: None,
            my_case: 1,
            my_first_s: 0.0,
            my_last_s: 1.0,
            my_first_u1: 0.0,
            my_last_u1: 0.0,
            my_first_u2: 0.0,
            my_last_u2: 0.0,
            my_length: 0.0,
            my_length1: 0.0,
            my_length2: 0.0,
            my_tol_len: the_tol,
            my_prev_s: 0.0,
            my_prev_u: 0.0,
            my_ui_1: None,
            my_si_1: None,
            my_ui_2: None,
            my_si_2: None,
        };
        this.init();
        this
    }

    /// OCCT Approx_CurvlinFunc(C2D, S, Tol) (cxx L139-152) — the case-2
    /// constructor.
    pub fn new_curve_on_surface(
        the_c2d: Curve2dHandle,
        the_s: SurfaceHandle,
        the_tol: f64,
    ) -> Self {
        let mut this = ApproxCurvlinFunc {
            my_c3d: None,
            my_c2d1: Some(the_c2d),
            my_c2d2: None,
            my_surf1: Some(the_s),
            my_surf2: None,
            my_case: 2,
            my_first_s: 0.0,
            my_last_s: 1.0,
            my_first_u1: 0.0,
            my_last_u1: 0.0,
            my_first_u2: 0.0,
            my_last_u2: 0.0,
            my_length: 0.0,
            my_length1: 0.0,
            my_length2: 0.0,
            my_tol_len: the_tol,
            my_prev_s: 0.0,
            my_prev_u: 0.0,
            my_ui_1: None,
            my_si_1: None,
            my_ui_2: None,
            my_si_2: None,
        };
        this.init();
        this
    }

    /// OCCT Approx_CurvlinFunc(C2D1, C2D2, S1, S2, Tol) (cxx L154-171) — the
    /// case-3 constructor.
    pub fn new_curve_on_2surfaces(
        the_c2d1: Curve2dHandle,
        the_c2d2: Curve2dHandle,
        the_s1: SurfaceHandle,
        the_s2: SurfaceHandle,
        the_tol: f64,
    ) -> Self {
        let mut this = ApproxCurvlinFunc {
            my_c3d: None,
            my_c2d1: Some(the_c2d1),
            my_c2d2: Some(the_c2d2),
            my_surf1: Some(the_s1),
            my_surf2: Some(the_s2),
            my_case: 3,
            my_first_s: 0.0,
            my_last_s: 1.0,
            my_first_u1: 0.0,
            my_last_u1: 0.0,
            my_first_u2: 0.0,
            my_last_u2: 0.0,
            my_length: 0.0,
            my_length1: 0.0,
            my_length2: 0.0,
            my_tol_len: the_tol,
            my_prev_s: 0.0,
            my_prev_u: 0.0,
            my_ui_1: None,
            my_si_1: None,
            my_ui_2: None,
            my_si_2: None,
        };
        this.init();
        this
    }

    /// OCCT Approx_CurvlinFunc::Init() (cxx L173-207).
    fn init(&mut self) {
        match self.my_case {
            1 => {
                // OCCT L180: Init(*myC3D, mySi_1, myUi_1).
                let a_c = self.my_c3d.as_ref().unwrap().clone();
                let (a_si, a_ui, a_first_u) = Self::init_curve(a_c.as_ref(), self.my_tol_len);
                // OCCT L261-262: the const_cast cache write of Init.
                self.my_prev_s = self.my_first_s;
                self.my_prev_u = a_first_u;
                self.my_si_1 = Some(a_si);
                self.my_ui_1 = Some(a_ui);
                // OCCT L181-183.
                self.my_first_u1 = self.my_c3d.as_ref().unwrap().first_parameter();
                self.my_last_u1 = self.my_c3d.as_ref().unwrap().last_parameter();
                self.my_first_u2 = 0.0;
                self.my_last_u2 = 0.0;
            }
            2 => {
                // OCCT L186-188: CurOnSur.Load(myC2D1); CurOnSur.Load(mySurf1);
                //               Init(CurOnSur, mySi_1, myUi_1).
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d1, a_surf1);
                let (a_si, a_ui, a_first_u) = Self::init_curve(&a_cur_on_sur, self.my_tol_len);
                // OCCT L261-262: the const_cast cache write of Init.
                self.my_prev_s = self.my_first_s;
                self.my_prev_u = a_first_u;
                self.my_si_1 = Some(a_si);
                self.my_ui_1 = Some(a_ui);
                // OCCT L189-191.
                self.my_first_u1 = Adaptor3dCurve::first_parameter(&a_cur_on_sur);
                self.my_last_u1 = Adaptor3dCurve::last_parameter(&a_cur_on_sur);
                self.my_first_u2 = 0.0;
                self.my_last_u2 = 0.0;
            }
            3 => {
                // OCCT L194-198.
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur1 =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d1, a_surf1);
                let (a_si1, a_ui1, a_first_u1) =
                    Self::init_curve(&a_cur_on_sur1, self.my_tol_len);
                // OCCT L261-262: the const_cast cache write of Init (the
                // first Init call, overwritten by the second below).
                self.my_prev_s = self.my_first_s;
                self.my_prev_u = a_first_u1;
                self.my_si_1 = Some(a_si1);
                self.my_ui_1 = Some(a_ui1);
                self.my_first_u1 = Adaptor3dCurve::first_parameter(&a_cur_on_sur1);
                self.my_last_u1 = Adaptor3dCurve::last_parameter(&a_cur_on_sur1);
                // OCCT L199-203.
                let a_c2d2 = self.my_c2d2.as_ref().unwrap().clone();
                let a_surf2 = self.my_surf2.as_ref().unwrap().clone();
                let a_cur_on_sur2 =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d2, a_surf2);
                let (a_si2, a_ui2, a_first_u2) =
                    Self::init_curve(&a_cur_on_sur2, self.my_tol_len);
                // OCCT L261-262 (the second Init call wins).
                self.my_prev_s = self.my_first_s;
                self.my_prev_u = a_first_u2;
                self.my_si_2 = Some(a_si2);
                self.my_ui_2 = Some(a_ui2);
                self.my_first_u2 = Adaptor3dCurve::first_parameter(&a_cur_on_sur2);
                self.my_last_u2 = Adaptor3dCurve::last_parameter(&a_cur_on_sur2);
            }
            _ => {}
        }

        // OCCT L206.
        self.length();
    }

    /// OCCT Approx_CurvlinFunc::Init(C, Si, Ui) (cxx L215-263) — the
    /// discretization arrays; the Si/Ui out-handles are the returned pair
    /// together with FirstU (the OCCT const_cast cache writes L261-262 are
    /// mirrored by the caller through (myPrevS, myPrevU)).
    fn init_curve(
        the_c: &dyn Adaptor3dCurve,
        the_tol_len: f64,
    ) -> (Vec<f64>, Vec<f64>, f64) {
        let a_first_u = the_c.first_parameter();
        let a_last_u = the_c.last_parameter();

        let a_nb_int: usize = 10;
        let a_nb_int_c3 = the_c.nb_intervals(GeomAbsShape::C3);
        // NCollection_Array1<double> Disc(1, NbIntC3 + 1).
        let mut a_disc = vec![0.0f64; a_nb_int_c3 + 1];

        if a_nb_int_c3 > 1 {
            let a_intervals = the_c.intervals(GeomAbsShape::C3);
            // OCCT C.Intervals(Disc, GeomAbs_C3) fills Disc(1..NbIntC3+1).
            for (a_i, a_v) in a_intervals.iter().enumerate() {
                if a_i < a_nb_int_c3 + 1 {
                    a_disc[a_i] = *a_v;
                }
            }
        } else {
            a_disc[0] = a_first_u;
            a_disc[1] = a_last_u;
        }

        let mut a_ui = vec![0.0f64; a_nb_int_c3 * a_nb_int + 1];
        let mut a_si = vec![0.0f64; a_nb_int_c3 * a_nb_int + 1];

        a_ui[0] = a_first_u;
        a_si[0] = 0.0;

        let mut a_i = 1usize;
        for a_j in 1..=a_nb_int_c3 {
            let a_step = (a_disc[a_j] - a_disc[a_j - 1]) / a_nb_int as f64;
            for _a_k in 1..=a_nb_int {
                a_ui[a_i] = a_ui[a_i - 1] + a_step;
                a_si[a_i] =
                    a_si[a_i - 1] + Self::length_curve(the_c, the_tol_len, a_ui[a_i - 1], a_ui[a_i]);
                a_i += 1;
            }
        }

        let a_len = a_si[a_si.len() - 1];
        for a_v in a_si.iter_mut() {
            *a_v /= a_len;
        }

        (a_si, a_ui, a_first_u)
    }

    /// OCCT Approx_CurvlinFunc::SetTol(Tol) (cxx L265-268).
    pub fn set_tol(&mut self, the_tol: f64) {
        self.my_tol_len = the_tol;
    }

    /// OCCT Approx_CurvlinFunc::FirstParameter() (cxx L270-273).
    pub fn first_parameter(&self) -> f64 {
        self.my_first_s
    }

    /// OCCT Approx_CurvlinFunc::LastParameter() (cxx L275-278).
    pub fn last_parameter(&self) -> f64 {
        self.my_last_s
    }

    /// OCCT Approx_CurvlinFunc::NbIntervals(S) (cxx L280-312).
    pub fn nb_intervals(&self, the_s: GeomAbsShape) -> usize {
        match self.my_case {
            1 => {
                // OCCT L287: return myC3D->NbIntervals(S).
                return self.my_c3d.as_ref().unwrap().nb_intervals(the_s);
            }
            2 => {
                // OCCT L289-291.
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d1, a_surf1);
                return a_cur_on_sur.nb_intervals(the_s);
            }
            3 => {
                // OCCT L293-307.
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur1 =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d1, a_surf1);
                let a_nb_int1 = a_cur_on_sur1.nb_intervals(the_s);
                let a_t1 = a_cur_on_sur1.intervals(the_s);
                let a_c2d2 = self.my_c2d2.as_ref().unwrap().clone();
                let a_surf2 = self.my_surf2.as_ref().unwrap().clone();
                let a_cur_on_sur2 =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d2, a_surf2);
                let _a_nb_int2 = a_cur_on_sur2.nb_intervals(the_s);
                let a_t2 = a_cur_on_sur2.intervals(the_s);

                let mut a_fusion: Vec<f64> = Vec::new();
                fuse_intervals(
                    &a_t1[..(a_nb_int1 + 1).min(a_t1.len())],
                    &a_t2,
                    &mut a_fusion,
                    // OCCT GeomLib::FuseIntervals(T1, T2, Fusion) — the 3-arg
                    // form (the GeomLib.hxx L209-210 defaults:
                    // Confusion = 1.0e-9, IsAdjustToFirstInterval = false).
                    1.0e-9,
                    false,
                );
                return a_fusion.len() - 1;
            }
            _ => {}
        }

        // OCCT L311: POP pour WNT.
        1
    }

    /// OCCT Approx_CurvlinFunc::Intervals(T, S) (cxx L314-355) — the T
    /// array values (the rcad Vec is the NCollection_Array1 result).
    pub fn intervals(&self, the_s: GeomAbsShape) -> Vec<f64> {
        let mut a_t: Vec<f64> = match self.my_case {
            1 => {
                // OCCT L322: myC3D->Intervals(T, S).
                self.my_c3d.as_ref().unwrap().intervals(the_s)
            }
            2 => {
                // OCCT L325-327.
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d1, a_surf1);
                a_cur_on_sur.intervals(the_s)
            }
            3 => {
                // OCCT L330-343.
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur1 =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d1, a_surf1);
                let a_nb_int1 = a_cur_on_sur1.nb_intervals(the_s);
                let a_t1 = a_cur_on_sur1.intervals(the_s);
                let a_c2d2 = self.my_c2d2.as_ref().unwrap().clone();
                let a_surf2 = self.my_surf2.as_ref().unwrap().clone();
                let a_cur_on_sur2 =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d2, a_surf2);
                let _a_nb_int2 = a_cur_on_sur2.nb_intervals(the_s);
                let a_t2 = a_cur_on_sur2.intervals(the_s);

                let mut a_fusion: Vec<f64> = Vec::new();
                // OCCT Approx_CurvlinFunc.cxx L338: GeomLib::FuseIntervals(T1,
                // T2, Fusion) — the GeomLib.hxx L206-210 defaults are
                // Confusion = 1.0e-9 and IsAdjustToFirstInterval = false
                // (not the PConfusion/true pair the nb_intervals site uses
                // with its own OCCT call shape).
                fuse_intervals(
                    &a_t1[..(a_nb_int1 + 1).min(a_t1.len())],
                    &a_t2,
                    &mut a_fusion,
                    1.0e-9,
                    false,
                );
                a_fusion
            }
            _ => Vec::new(),
        };

        // OCCT L351-354.
        for a_v in a_t.iter_mut() {
            *a_v = self.get_s_parameter(*a_v);
        }
        a_t
    }

    /// OCCT Approx_CurvlinFunc::Trim(First, Last, Tol) (cxx L357-413).
    pub fn trim(&mut self, the_first: f64, the_last: f64, the_tol: f64) {
        if the_first < 0.0 || the_last > 1.0 {
            // throw Standard_OutOfRange("Approx_CurvlinFunc::Trim").
            panic!("Standard_OutOfRange: Approx_CurvlinFunc::Trim");
        }
        if (the_last - the_first) < the_tol {
            return;
        }

        match self.my_case {
            1 => {
                // OCCT L375: myC3D = myC3D->Trim(myFirstU1, myLastU1, Tol).
                let a_c3d = self.my_c3d.take().unwrap();
                let a_c3d = a_c3d.trim(self.my_first_u1, self.my_last_u1, the_tol);
                // OCCT L376-377.
                let a_first_u = self.get_u_parameter(a_c3d.as_ref(), the_first, 1);
                let a_last_u = self.get_u_parameter(a_c3d.as_ref(), the_last, 1);
                // OCCT L378: myC3D = myC3D->Trim(FirstU, LastU, Tol).
                self.my_c3d = Some(a_c3d.trim(a_first_u, a_last_u, the_tol));
            }
            3 => {
                // OCCT L381-393 (case-3: the second curve is trimmed first).
                {
                    let a_c2d2 = self.my_c2d2.as_ref().unwrap().clone();
                    let a_surf2 = self.my_surf2.as_ref().unwrap().clone();
                    let a_cur_on_sur =
                        rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(
                            a_c2d2.clone(),
                            a_surf2.clone(),
                        );
                    // OCCT L383: HCurOnSur = down_cast<...>(CurOnSur.Trim(
                    //               myFirstU2, myLastU2, Tol)); the rcad
                    // trim encoding rebuilds the (curve, surface) pair
                    // (adaptor.rs L1609-1612) — the surface is unchanged.
                    let a_c2d2_trim = a_c2d2.trim(self.my_first_u2, self.my_last_u2, the_tol);
                    self.my_c2d2 = Some(a_c2d2_trim);
                    self.my_surf2 = Some(a_surf2);
                    let a_cur_on_sur =
                        rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(
                            self.my_c2d2.as_ref().unwrap().clone(),
                            self.my_surf2.as_ref().unwrap().clone(),
                        );
                    // OCCT L389-393.
                    let a_first_u = self.get_u_parameter(&a_cur_on_sur, the_first, 1);
                    let a_last_u = self.get_u_parameter(&a_cur_on_sur, the_last, 1);
                    let a_c2d2 = self.my_c2d2.as_ref().unwrap().clone();
                    self.my_c2d2 = Some(a_c2d2.trim(a_first_u, a_last_u, the_tol));
                }
                // OCCT L395 [[fallthrough]] into case 2.
                {
                    let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                    let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                    let a_cur_on_sur =
                        rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(
                            a_c2d1.clone(),
                            a_surf1.clone(),
                        );
                    let a_c2d1_trim = a_c2d1.trim(self.my_first_u1, self.my_last_u1, the_tol);
                    self.my_c2d1 = Some(a_c2d1_trim);
                    self.my_surf1 = Some(a_surf1);
                    let a_cur_on_sur =
                        rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(
                            self.my_c2d1.as_ref().unwrap().clone(),
                            self.my_surf1.as_ref().unwrap().clone(),
                        );
                    let a_first_u = self.get_u_parameter(&a_cur_on_sur, the_first, 1);
                    let a_last_u = self.get_u_parameter(&a_cur_on_sur, the_last, 1);
                    let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                    self.my_c2d1 = Some(a_c2d1.trim(a_first_u, a_last_u, the_tol));
                }
            }
            2 => {
                // OCCT L397-409.
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(
                        a_c2d1.clone(),
                        a_surf1.clone(),
                    );
                let a_c2d1_trim = a_c2d1.trim(self.my_first_u1, self.my_last_u1, the_tol);
                self.my_c2d1 = Some(a_c2d1_trim);
                self.my_surf1 = Some(a_surf1);
                let a_cur_on_sur = rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(
                    self.my_c2d1.as_ref().unwrap().clone(),
                    self.my_surf1.as_ref().unwrap().clone(),
                );
                let a_first_u = self.get_u_parameter(&a_cur_on_sur, the_first, 1);
                let a_last_u = self.get_u_parameter(&a_cur_on_sur, the_last, 1);
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                self.my_c2d1 = Some(a_c2d1.trim(a_first_u, a_last_u, the_tol));
            }
            _ => {}
        }
        // OCCT L411-412.
        self.my_first_s = the_first;
        self.my_last_s = the_last;
    }

    /// OCCT Approx_CurvlinFunc::Length() (cxx L415-449).
    pub fn length(&mut self) {
        match self.my_case {
            1 => {
                // OCCT L423-426.
                let a_first_u = self.my_c3d.as_ref().unwrap().first_parameter();
                let a_last_u = self.my_c3d.as_ref().unwrap().last_parameter();
                let a_c = self.my_c3d.as_ref().unwrap().clone();
                self.my_length = Self::length_curve(a_c.as_ref(), self.my_tol_len, a_first_u, a_last_u);
                self.my_length1 = 0.0;
                self.my_length2 = 0.0;
            }
            2 => {
                // OCCT L429-434.
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d1, a_surf1);
                let a_first_u = Adaptor3dCurve::first_parameter(&a_cur_on_sur);
                let a_last_u = Adaptor3dCurve::last_parameter(&a_cur_on_sur);
                self.my_length =
                    Self::length_curve(&a_cur_on_sur, self.my_tol_len, a_first_u, a_last_u);
                self.my_length1 = 0.0;
                self.my_length2 = 0.0;
            }
            3 => {
                // OCCT L437-447.
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur1 =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d1, a_surf1);
                let a_first_u = Adaptor3dCurve::first_parameter(&a_cur_on_sur1);
                let a_last_u = Adaptor3dCurve::last_parameter(&a_cur_on_sur1);
                self.my_length1 =
                    Self::length_curve(&a_cur_on_sur1, self.my_tol_len, a_first_u, a_last_u);
                let a_c2d2 = self.my_c2d2.as_ref().unwrap().clone();
                let a_surf2 = self.my_surf2.as_ref().unwrap().clone();
                let a_cur_on_sur2 =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d2, a_surf2);
                let a_first_u = Adaptor3dCurve::first_parameter(&a_cur_on_sur2);
                let a_last_u = Adaptor3dCurve::last_parameter(&a_cur_on_sur2);
                self.my_length2 =
                    Self::length_curve(&a_cur_on_sur2, self.my_tol_len, a_first_u, a_last_u);
                self.my_length = (self.my_length1 + self.my_length2) / 2.0;
            }
            _ => {}
        }
    }

    /// OCCT Approx_CurvlinFunc::Length(C, FirstU, LastU) (cxx L451-457).
    pub fn length_curve(
        the_c: &dyn Adaptor3dCurve,
        the_tol_len: f64,
        the_first_u: f64,
        the_last_u: f64,
    ) -> f64 {
        // OCCT L455: GCPnts_AbscissaPoint::Length(C, FirstU, LastU, myTolLen).
        gcpnts_length_3d_range_tol(&AdaptorAsGCPnts(the_c), the_first_u, the_last_u, the_tol_len)
    }

    /// OCCT Approx_CurvlinFunc::GetLength() (cxx L459-462).
    pub fn get_length(&self) -> f64 {
        self.my_length
    }

    /// OCCT Approx_CurvlinFunc::GetSParameter(U) (cxx L464-489).
    pub fn get_s_parameter(&self, the_u: f64) -> f64 {
        let mut a_s;
        match self.my_case {
            1 => {
                // OCCT L472.
                let a_c = self.my_c3d.as_ref().unwrap().clone();
                a_s = Self::get_s_parameter_curve(
                    a_c.as_ref(),
                    self.my_tol_len,
                    self.my_first_s,
                    the_u,
                    self.my_length,
                );
            }
            2 => {
                // OCCT L475-477.
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d1, a_surf1);
                a_s = Self::get_s_parameter_curve(
                    &a_cur_on_sur,
                    self.my_tol_len,
                    self.my_first_s,
                    the_u,
                    self.my_length,
                );
            }
            3 => {
                // OCCT L480-486.
                let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
                let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
                let a_cur_on_sur1 =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d1, a_surf1);
                let a_s1 = Self::get_s_parameter_curve(
                    &a_cur_on_sur1,
                    self.my_tol_len,
                    self.my_first_s,
                    the_u,
                    self.my_length1,
                );
                let a_c2d2 = self.my_c2d2.as_ref().unwrap().clone();
                let a_surf2 = self.my_surf2.as_ref().unwrap().clone();
                let a_cur_on_sur2 =
                    rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(a_c2d2, a_surf2);
                let a_s2 = Self::get_s_parameter_curve(
                    &a_cur_on_sur2,
                    self.my_tol_len,
                    self.my_first_s,
                    the_u,
                    self.my_length2,
                );
                a_s = (a_s1 + a_s2) / 2.0;
            }
            _ => {
                a_s = 0.0;
            }
        }
        a_s
    }

    /// OCCT Approx_CurvlinFunc::GetUParameter(C, S, NumberOfCurve)
    /// (cxx L491-572).
    pub fn get_u_parameter(
        &mut self,
        the_c: &dyn Adaptor3dCurve,
        the_s: f64,
        the_number_of_curve: i32,
    ) -> f64 {
        if the_s < 0.0 || the_s > 1.0 {
            // throw Standard_ConstructionError("Approx_CurvlinFunc::GetUParameter").
            panic!("Standard_ConstructionError: Approx_CurvlinFunc::GetUParameter");
        }

        let (the_init_u_array, the_init_s_array, the_length) = if the_number_of_curve == 1 {
            (
                self.my_ui_1.as_ref().unwrap().clone(),
                self.my_si_1.as_ref().unwrap().clone(),
                if self.my_case == 3 {
                    self.my_length1
                } else {
                    self.my_length
                },
            )
        } else {
            (
                self.my_ui_2.as_ref().unwrap().clone(),
                self.my_si_2.as_ref().unwrap().clone(),
                self.my_length2,
            )
        };

        let a_nb_int = the_init_u_array.len() as i32 - 1;

        let mut a_n_interval;
        if the_s == 1.0 {
            a_n_interval = a_nb_int - 1;
        } else {
            let mut a_i = 0i32;
            while a_i < a_nb_int {
                // OCCT L536: the Value(i) reads are 0-based on the (0, N)
                // arrays — direct indexing, not the 1-based `at` helper.
                if the_init_s_array[a_i as usize] <= the_s
                    && the_s < the_init_s_array[(a_i + 1) as usize]
                {
                    break;
                }
                a_i += 1;
            }
            a_n_interval = a_i;
        }
        if the_s == the_init_s_array[a_n_interval as usize] {
            return the_init_u_array[a_n_interval as usize];
        }
        if the_s == the_init_s_array[(a_n_interval + 1) as usize] {
            return the_init_u_array[(a_n_interval + 1) as usize];
        }

        let a_base = the_init_u_array[a_n_interval as usize];
        let a_delta_s = (the_s - the_init_s_array[a_n_interval as usize]) * the_length;

        // to find an initial point
        let mut a_xi = [0.0f64; 4];
        let mut a_yi = [0.0f64; 4];
        findfourpoints(
            the_s,
            a_n_interval,
            &the_init_s_array,
            &the_init_u_array,
            self.my_prev_s,
            self.my_prev_u,
            &mut a_xi,
            &mut a_yi,
        );
        let a_uguess = cubic(the_s, &a_xi, &a_yi);

        // OCCT L560: U = GCPnts_AbscissaPoint(C, deltaS, base, UGuess,
        //               myTolLen).Parameter().
        let the_u = gcpnts_abscissa_point_adv(the_c, a_delta_s, a_base, a_uguess, self.my_tol_len);

        // OCCT L563-564: the const_cast cache write.
        self.my_prev_s = the_s;
        self.my_prev_u = the_u;

        the_u
    }

    /// OCCT Approx_CurvlinFunc::GetSParameter(C, U, Len) (cxx L574-581).
    fn get_s_parameter_curve(
        the_c: &dyn Adaptor3dCurve,
        the_tol_len: f64,
        the_first_s: f64,
        the_u: f64,
        the_len: f64,
    ) -> f64 {
        let a_origin = the_c.first_parameter();
        let a_s = the_first_s + Self::length_curve(the_c, the_tol_len, a_origin, the_u) / the_len;
        a_s
    }

    /// OCCT Approx_CurvlinFunc::EvalCase1(S, Order, Result) (cxx L583-637).
    pub fn eval_case1(&mut self, the_s: f64, the_order: i32, the_result: &mut [f64]) -> bool {
        if self.my_case != 1 {
            panic!("Standard_ConstructionError: Approx_CurvlinFunc::EvalCase1");
        }

        let a_c3d = self.my_c3d.as_ref().unwrap().clone();
        let the_u = self.get_u_parameter(a_c3d.as_ref(), the_s, 1);

        match the_order {
            0 => {
                // OCCT L602: myC3D->D0(U, C).
                let a_c = a_c3d.value(the_u);

                the_result[0] = a_c.x;
                the_result[1] = a_c.y;
                the_result[2] = a_c.z;
            }
            1 => {
                // OCCT L610-617.
                let (a_c, a_dc_du) = a_c3d.d1(the_u);
                let a_mag = a_dc_du.length();
                let a_du_ds = self.my_length / a_mag;
                let a_dc_ds = a_dc_du * a_du_ds;

                the_result[0] = a_dc_ds.x;
                the_result[1] = a_dc_ds.y;
                the_result[2] = a_dc_ds.z;
                let _ = a_c;
            }
            2 => {
                // OCCT L621-630.
                let (_a_c, a_dc_du, a_d2c_du2) = a_c3d.d2(the_u);
                let a_mag = a_dc_du.length();
                let a_du_ds = self.my_length / a_mag;
                let a_d2u_ds2 = -self.my_length * a_dc_du.dot(a_d2c_du2) * a_du_ds / (a_mag * a_mag * a_mag);
                let a_d2c_ds2 = a_d2c_du2 * (a_du_ds * a_du_ds) + a_dc_du * a_d2u_ds2;

                the_result[0] = a_d2c_ds2.x;
                the_result[1] = a_d2c_ds2.y;
                the_result[2] = a_d2c_ds2.z;
            }
            _ => {
                // OCCT L632-634.
                the_result[0] = 0.0;
                the_result[1] = 0.0;
                the_result[2] = 0.0;
                return false;
            }
        }
        true
    }

    /// OCCT Approx_CurvlinFunc::EvalCase2(S, Order, Result) (cxx L639-653).
    pub fn eval_case2(&mut self, the_s: f64, the_order: i32, the_result: &mut [f64]) -> bool {
        if self.my_case != 2 {
            panic!("Standard_ConstructionError: Approx_CurvlinFunc::EvalCase2");
        }

        // OCCT L650: Done = EvalCurOnSur(S, Order, Result, 1).
        let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
        let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
        self.eval_cur_on_sur(the_s, the_order, the_result, 1, &a_c2d1, &a_surf1)
    }

    /// OCCT Approx_CurvlinFunc::EvalCase3(S, Order, Result) (cxx L655-680).
    pub fn eval_case3(&mut self, the_s: f64, the_order: i32, the_result: &mut [f64]) -> bool {
        if self.my_case != 3 {
            panic!("Standard_ConstructionError: Approx_CurvlinFunc::EvalCase3");
        }

        let mut a_tmp_res1 = [0.0f64; 5];
        let mut a_tmp_res2 = [0.0f64; 5];

        let a_c2d1 = self.my_c2d1.as_ref().unwrap().clone();
        let a_surf1 = self.my_surf1.as_ref().unwrap().clone();
        let mut a_done = self.eval_cur_on_sur(the_s, the_order, &mut a_tmp_res1, 1, &a_c2d1, &a_surf1);

        let a_c2d2 = self.my_c2d2.as_ref().unwrap().clone();
        let a_surf2 = self.my_surf2.as_ref().unwrap().clone();
        a_done = self.eval_cur_on_sur(the_s, the_order, &mut a_tmp_res2, 2, &a_c2d2, &a_surf2) && a_done;

        // OCCT L671-677.
        the_result[0] = a_tmp_res1[0];
        the_result[1] = a_tmp_res1[1];
        the_result[2] = a_tmp_res2[0];
        the_result[3] = a_tmp_res2[1];
        the_result[4] = 0.5 * (a_tmp_res1[2] + a_tmp_res2[2]);
        the_result[5] = 0.5 * (a_tmp_res1[3] + a_tmp_res2[3]);
        the_result[6] = 0.5 * (a_tmp_res1[4] + a_tmp_res2[4]);

        a_done
    }

    /// OCCT Approx_CurvlinFunc::EvalCurOnSur(S, Order, Result, NumberOfCurve)
    /// (cxx L682-793) — the 2d/3d evaluation on one of the two surfaces.
    fn eval_cur_on_sur(
        &mut self,
        the_s: f64,
        the_order: i32,
        the_result: &mut [f64],
        the_number_of_curve: i32,
        the_c2d: &Curve2dHandle,
        the_surf: &SurfaceHandle,
    ) -> bool {
        // OCCT L691-717: the curve/surface pair and the length of the
        // selected curve.
        let the_length = if the_number_of_curve == 1 {
            if self.my_case == 3 {
                self.my_length1
            } else {
                self.my_length
            }
        } else if the_number_of_curve == 2 {
            self.my_length2
        } else {
            panic!("Standard_ConstructionError: Approx_CurvlinFunc::EvalCurOnSur");
        };
        let a_cur_on_sur = rcad_kernel::base::proj_lib::adaptor::CurveOnSurface::new(
            the_c2d.clone(),
            the_surf.clone(),
        );
        let the_u = self.get_u_parameter(&a_cur_on_sur, the_s, the_number_of_curve);

        match the_order {
            0 => {
                // OCCT L728-735.
                let a_c2d = the_c2d.d0(the_u);
                let a_c = the_surf.value(a_c2d.x, a_c2d.y);

                the_result[0] = a_c2d.x;
                the_result[1] = a_c2d.y;
                the_result[2] = a_c.x;
                the_result[3] = a_c.y;
                the_result[4] = a_c.z;
            }
            1 => {
                // OCCT L739-755.
                let (a_c2d, a_dc2d_du) = the_c2d.d1(the_u);
                let a_dv_du = a_dc2d_du.x;
                let a_dw_du = a_dc2d_du.y;
                let (a_c, a_ds_dv, a_ds_dw) = the_surf.d1(a_c2d.x, a_c2d.y);
                let a_dc_du = a_ds_dv * a_dv_du + a_ds_dw * a_dw_du;
                let a_mag = a_dc_du.length();
                let a_du_ds = the_length / a_mag;

                let a_dv_ds = a_dv_du * a_du_ds;
                let a_dw_ds = a_dw_du * a_du_ds;
                let a_dc_ds = a_dc_du * a_du_ds;

                the_result[0] = a_dv_ds;
                the_result[1] = a_dw_ds;
                the_result[2] = a_dc_ds.x;
                the_result[3] = a_dc_ds.y;
                the_result[4] = a_dc_ds.z;
                let _ = a_c;
            }
            2 => {
                // OCCT L759-785.
                let (a_c2d, a_dc2d_du, a_d2c2d_du2) = the_c2d.d2(the_u);
                let a_dv_du = a_dc2d_du.x;
                let a_dw_du = a_dc2d_du.y;
                let a_d2v_du2 = a_d2c2d_du2.x;
                let a_d2w_du2 = a_d2c2d_du2.y;
                let (_a_c, a_ds_dv, a_ds_dw, a_d2s_dv2, a_d2s_dw2, a_d2s_dvdw) =
                    the_surf.d2(a_c2d.x, a_c2d.y);
                let a_dc_du = a_ds_dv * a_dv_du + a_ds_dw * a_dw_du;
                let a_d2c_du2 = (a_d2s_dv2 * a_dv_du + a_d2s_dvdw * a_dw_du) * a_dv_du
                    + a_ds_dv * a_d2v_du2
                    + (a_d2s_dvdw * a_dv_du + a_d2s_dw2 * a_dw_du) * a_dw_du
                    + a_ds_dw * a_d2w_du2;
                let a_mag = a_dc_du.length();
                let a_du_ds = the_length / a_mag;
                // OCCT L770: d2U_dS2 = -Length * C'.C'' * dU_dS / Mag^3 —
                // the value consumed by the d2V_dS2 / d2W_dS2 lines below.
                let a_d2u_ds2 =
                    -the_length * a_dc_du.dot(a_d2c_du2) * a_du_ds / (a_mag * a_mag * a_mag);

                // OCCT L772-775.
                let a_dv_ds = a_dv_du * a_du_ds;
                let a_dw_ds = a_dw_du * a_du_ds;
                let a_d2v_ds2 = a_d2v_du2 * (a_du_ds * a_du_ds) + a_dv_du * a_d2u_ds2;
                let a_d2w_ds2 = a_d2w_du2 * (a_du_ds * a_du_ds) + a_dw_du * a_d2u_ds2;

                // OCCT L777: d2U_dS2 = -C'.C'' * dU_dS / Mag^2 — the OCCT
                // reassignment (the value is not read afterwards; kept as
                // the recorded dead write).
                let _a_d2u_ds2 = -a_dc_du.dot(a_d2c_du2) * a_du_ds / (a_mag * a_mag);
                // OCCT L778-779.
                let a_d2c_ds2 = (a_d2s_dv2 * a_dv_ds + a_d2s_dvdw * a_dw_ds) * a_dv_ds
                    + a_ds_dv * a_d2v_ds2
                    + (a_d2s_dw2 * a_dw_ds + a_d2s_dvdw * a_dv_ds) * a_dw_ds
                    + a_ds_dw * a_d2w_ds2;

                the_result[0] = a_d2v_ds2;
                the_result[1] = a_d2w_ds2;
                the_result[2] = a_d2c_ds2.x;
                the_result[3] = a_d2c_ds2.y;
                the_result[4] = a_d2c_ds2.z;
            }
            _ => {
                // OCCT L788-790.
                the_result[0] = 0.0;
                the_result[1] = 0.0;
                the_result[2] = 0.0;
                the_result[3] = 0.0;
                the_result[4] = 0.0;
                return false;
            }
        }
        true
    }
}

// =========================================================================
// Tests (hand-derived).
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::base::proj_lib::GeomCurveAdaptor;
    use rcad_kernel::geom::{Curve3, Line3};

    /// Hand derivation: the reference curve is the X segment u in [0, 2]
    /// (|C'(u)| = 1, total length 2).  The curvilinear parameter is the
    /// normalized abscissa: S(u) = u / 2 and U(S) = 2 * S.
    fn segment_curve() -> CurveHandle {
        Arc::new(GeomCurveAdaptor::with_range(
            Curve3::Line(Line3::new(glam::DVec3::ZERO, glam::DVec3::X)),
            0.0,
            2.0,
        ))
    }

    #[test]
    fn case1_segment_s_parameter_and_length() {
        let mut a_f = ApproxCurvlinFunc::new_curve(segment_curve(), 1.0e-6);
        // The parameter window is the normalized abscissa [0, 1].
        assert_eq!(a_f.first_parameter(), 0.0);
        assert_eq!(a_f.last_parameter(), 1.0);
        // The stored Si array is normalized: Si(Upper) == 1.
        let a_si = a_f.my_si_1.as_ref().unwrap();
        assert!((a_si[a_si.len() - 1] - 1.0).abs() < 1.0e-9);
        // Length(C, 0, 2) = 2 for the unit-speed segment.
        assert!((a_f.get_length() - 2.0).abs() < 1.0e-9);
        // GetSParameter(1.5) = myFirstS + Length(0, 1.5) / 2 = 0.75.
        let a_s = a_f.get_s_parameter(1.5);
        assert!((a_s - 0.75).abs() < 1.0e-9, "S(1.5) = {}", a_s);
    }

    #[test]
    fn case1_segment_get_u_parameter_inverts_s() {
        let mut a_f = ApproxCurvlinFunc::new_curve(segment_curve(), 1.0e-6);
        let a_c = segment_curve();
        // U(0.75) must land back on u = 1.5 (Newton on Integral(0, u, 1)).
        let a_u = a_f.get_u_parameter(a_c.as_ref(), 0.75, 1);
        assert!((a_u - 1.5).abs() < 1.0e-6, "U(0.75) = {}", a_u);
        // The knots: U(0) = 0, U(1) = 2.
        let a_u0 = a_f.get_u_parameter(a_c.as_ref(), 0.0, 1);
        let a_u1 = a_f.get_u_parameter(a_c.as_ref(), 1.0, 1);
        assert!((a_u0 - 0.0).abs() < 1.0e-12);
        assert!((a_u1 - 2.0).abs() < 1.0e-9);
    }

    #[test]
    fn case1_segment_eval_derivatives() {
        let mut a_f = ApproxCurvlinFunc::new_curve(segment_curve(), 1.0e-6);
        // Order 0 at S = 0.5: C(U(0.5)) = C(1.0) = (1, 0, 0).
        let mut a_r0 = [0.0f64; 3];
        assert!(a_f.eval_case1(0.5, 0, &mut a_r0));
        assert!((a_r0[0] - 1.0).abs() < 1.0e-9 && a_r0[1].abs() < 1.0e-9 && a_r0[2].abs() < 1.0e-9);
        // Order 1: dC/dS = C'(u) * dU/dS with dU/dS = Length / |C'| = 2.
        let mut a_r1 = [0.0f64; 3];
        assert!(a_f.eval_case1(0.5, 1, &mut a_r1));
        assert!(
            (a_r1[0] - 2.0).abs() < 1.0e-9 && a_r1[1].abs() < 1.0e-9 && a_r1[2].abs() < 1.0e-9,
            "dC/dS = {:?}",
            a_r1
        );
        // Order 2: C'' = 0 on the segment, so d2C/dS2 = 0.
        let mut a_r2 = [0.0f64; 3];
        assert!(a_f.eval_case1(0.5, 2, &mut a_r2));
        assert!(a_r2.iter().all(|a_v| a_v.abs() < 1.0e-9));
        // Order 3 is not delivered (the OCCT default arm answers false).
        let mut a_r3 = [9.0f64; 3];
        assert!(!a_f.eval_case1(0.5, 3, &mut a_r3));
        assert!(a_r3.iter().all(|a_v| *a_v == 0.0));
    }

    #[test]
    fn cubic_and_findfourpoints_hand_values() {
        // cubic over the identity points (0,0) (1,1) (2,2) (3,3) at X = 2.5.
        let a_xi = [0.0f64, 1.0, 2.0, 3.0];
        let a_yi = [0.0f64, 1.0, 2.0, 3.0];
        assert!((cubic(2.5, &a_xi, &a_yi) - 2.5).abs() < 1.0e-12);
        // findfourpoints clamps NInterval into [1, NbInt - 2] and picks the
        // four knots (S, U) starting at NInterval - 1.
        let a_si = vec![0.0f64, 0.25, 0.5, 0.75, 1.0];
        let a_ui = vec![0.0f64, 0.5, 1.0, 1.5, 2.0];
        let mut a_xi2 = [0.0f64; 4];
        let mut a_yi2 = [0.0f64; 4];
        findfourpoints(0.6, 9, &a_si, &a_ui, 5.0, 5.0, &mut a_xi2, &mut a_yi2);
        // NInterval = 9 clamps to NbInt - 2 = 2 -> knots Si(1..4) 0-based;
        // prevS = 5 falls outside every (Xi[i], Xi[i+1]) — no insert.
        assert_eq!(a_xi2, [0.25, 0.5, 0.75, 1.0]);
        assert_eq!(a_yi2, [0.5, 1.0, 1.5, 2.0]);
    }
}
