//! OCCT GccEnt + Geom2dGcc — 2D constraint construction (tangent line to a
//! curve and a point / to two curves).
//!
//! 1:1 translations:
//!   - GccEnt_Position.hxx — the position qualifier enum.
//!   - Geom2dGcc_QualifiedCurve.hxx/.cxx + Geom2dGcc_QCurve.hxx/.cxx.
//!   - Geom2dGcc_CurveTool.hxx/.cxx — the 2D adaptor tool (FirstParameter/
//!     LastParameter/NbSamples/EpsX/D1/D2/D3).
//!   - Geom2dGcc_FunctionTanCuPnt.cxx — F(u) = P1P2×T/(|P1P2|·|T|) tangency
//!     function (curve vs point).
//!   - Geom2dGcc_FunctionTanCuCu.cxx — the two-curve tangency function set.
//!   - Geom2dGcc_Lin2d2TanIter.cxx (whole) — the iterative tangent-line solver
//!     driven by math_FunctionRoot / math_FunctionSetRoot.
//!   - Geom2dGcc_Lin2d2Tan.cxx (whole) — the dispatcher (analytic circle case
//!     via GccAna_Lin2d2Tan for (circle, point/circle) pairs — ported for the
//!     point-point/circle cases the tests hit; general curves via the sampling
//!     loop + Lin2d2TanIter + the Add dedup) plus the result accessors.

use glam::DVec2;
use rcad_kernel::core::precision::parametric_default;
use rcad_kernel::geom::{Curve2d, Curve2dEval};
use rcad_kernel::math::function_set_root::{FunctionSetRoot, FunctionSetWithDerivatives};

/// OCCT GccEnt_Position (GccEnt_Position.hxx).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GccEntPosition {
    Unqualified,
    Enclosing,
    Enclosed,
    Outside,
    NoQualifier,
}

/// OCCT Geom2dGcc_QualifiedCurve — a 2D curve with a position qualifier.
#[derive(Debug, Clone)]
pub struct QualifiedCurve {
    qualified: Curve2d,
    qualifier: GccEntPosition,
}

impl QualifiedCurve {
    /// OCCT Geom2dGcc_QualifiedCurve(Curve, Qualifier).
    pub fn new(curve: Curve2d, qualifier: GccEntPosition) -> Self {
        QualifiedCurve {
            qualified: curve,
            qualifier,
        }
    }

    /// OCCT Qualified().
    pub fn qualified(&self) -> &Curve2d {
        &self.qualified
    }

    /// OCCT Qualifier().
    pub fn qualifier(&self) -> GccEntPosition {
        self.qualifier
    }

    /// OCCT IsUnqualified().
    pub fn is_unqualified(&self) -> bool {
        self.qualifier == GccEntPosition::Unqualified
    }
    /// OCCT IsEnclosing().
    pub fn is_enclosing(&self) -> bool {
        self.qualifier == GccEntPosition::Enclosing
    }
    /// OCCT IsEnclosed().
    pub fn is_enclosed(&self) -> bool {
        self.qualifier == GccEntPosition::Enclosed
    }
    /// OCCT IsOutside().
    pub fn is_outside(&self) -> bool {
        self.qualifier == GccEntPosition::Outside
    }
}

/// OCCT Geom2dGcc_QCurve — the adaptor + qualifier used by the iterative
/// solver.
#[derive(Debug, Clone)]
pub struct QCurve {
    qualified: Curve2d,
    qualifier: GccEntPosition,
}

impl QCurve {
    /// OCCT Geom2dGcc_QCurve(Curve, Qualifier).
    pub fn new(curve: Curve2d, qualifier: GccEntPosition) -> Self {
        QCurve {
            qualified: curve,
            qualifier,
        }
    }

    /// OCCT Qualified().
    pub fn qualified(&self) -> &Curve2d {
        &self.qualified
    }

    /// OCCT Qualifier().
    pub fn qualifier(&self) -> GccEntPosition {
        self.qualifier
    }

    /// OCCT IsUnqualified() etc.
    pub fn is_unqualified(&self) -> bool {
        self.qualifier == GccEntPosition::Unqualified
    }
    pub fn is_enclosing(&self) -> bool {
        self.qualifier == GccEntPosition::Enclosing
    }
    pub fn is_enclosed(&self) -> bool {
        self.qualifier == GccEntPosition::Enclosed
    }
    pub fn is_outside(&self) -> bool {
        self.qualifier == GccEntPosition::Outside
    }
}

/// OCCT Geom2dGcc_CurveTool (Geom2dGcc_CurveTool.cxx) — the adaptor tool.
pub mod curve_tool {
    use super::*;

    /// OCCT Geom2dGcc_CurveTool::FirstParameter.
    pub fn first_parameter(c: &Curve2d) -> f64 {
        c.default_domain()[0]
    }

    /// OCCT Geom2dGcc_CurveTool::LastParameter.
    pub fn last_parameter(c: &Curve2d) -> f64 {
        c.default_domain()[1]
    }

    /// OCCT Geom2dGcc_CurveTool::NbSamples — 20.
    pub fn nb_samples(_c: &Curve2d) -> usize {
        20
    }

    /// OCCT Geom2dGcc_CurveTool::EpsX — the adaptor Resolution(Tol):
    /// line → Tol, circle → 2·asin(Tol/2R), other conics →
    /// Precision::Parametric(Tol) = Tol·0.01.
    pub fn eps_x(c: &Curve2d, tol: f64) -> f64 {
        match c {
            Curve2d::Line(_) => tol,
            Curve2d::Circle(cir) => {
                let r = cir.radius;
                if r > tol / 2.0 {
                    2.0 * (tol / (2.0 * r)).asin()
                } else {
                    std::f64::consts::TAU
                }
            }
            _ => parametric_default(tol),
        }
    }

    /// OCCT Geom2dGcc_CurveTool::Value.
    pub fn value(c: &Curve2d, u: f64) -> DVec2 {
        c.point_at(u)
    }

    /// OCCT Geom2dGcc_CurveTool::D1(U, P, T).
    pub fn d1(c: &Curve2d, u: f64) -> (DVec2, DVec2) {
        (c.point_at(u), c.derivative_at(u))
    }

    /// OCCT Geom2dGcc_CurveTool::D2(U, P, T, N).
    pub fn d2(c: &Curve2d, u: f64) -> (DVec2, DVec2, DVec2) {
        (c.point_at(u), c.derivative_at(u), c.derivative2_at(u))
    }
}

/// OCCT Geom2dGcc_FunctionTanCuPnt (Geom2dGcc_FunctionTanCuPnt.cxx) —
/// F(u) = (Point→P(u))×T(u) / (|T(u)|·|Point→P(u)|).
#[derive(Debug, Clone)]
pub struct FunctionTanCuPnt {
    curve: Curve2d,
    point: DVec2,
}

impl FunctionTanCuPnt {
    /// OCCT Geom2dGcc_FunctionTanCuPnt(C, Point).
    pub fn new(curve: Curve2d, point: DVec2) -> Self {
        FunctionTanCuPnt { curve, point }
    }

    /// OCCT Value(X, Fval).
    pub fn value(&self, x: f64) -> f64 {
        let (point, vect) = curve_tool::d1(&self.curve, x);
        let norme_d1 = vect.length();
        let the_direction = point - self.point;
        let norme_dir = the_direction.length();
        cross2d(the_direction, vect) / (norme_d1 * norme_dir)
    }

    /// OCCT Derivative(X, Deriv).
    pub fn derivative(&self, x: f64) -> f64 {
        let (point, vec1, vec2) = curve_tool::d2(&self.curve, x);
        let the_direction = point - self.point;
        let norme_d1 = vec1.length();
        let norme_dir = the_direction.length();
        cross2d(the_direction, vec2) / (norme_d1 * norme_dir)
            - (cross2d(the_direction, vec1) / (norme_d1 * norme_dir))
                * (vec1.dot(vec2) / (norme_d1 * norme_d1)
                    + vec1.dot(the_direction) / (norme_dir * norme_dir))
    }
}

/// OCCT Geom2dGcc_FunctionTanCuCu (Geom2dGcc_FunctionTanCuCu.cxx) — the
/// two-curve tangency function set:
///   F1 = P1P2×T1 / (|T1|·|P1P2|²),  F2 = T1×T2 / (|T1|·|T2|).
#[derive(Debug, Clone)]
pub struct FunctionTanCuCu {
    curve1: Curve2d,
    curve2: Curve2d,
}

impl FunctionTanCuCu {
    /// OCCT Geom2dGcc_FunctionTanCuCu(C1, C2).
    pub fn new(curve1: Curve2d, curve2: Curve2d) -> Self {
        FunctionTanCuCu { curve1, curve2 }
    }

    /// OCCT InitDerivative(X, Point1, Point2, Tan1, Tan2, D21, D22).
    fn init_derivative(&self, x: &[f64]) -> (DVec2, DVec2, DVec2, DVec2, DVec2, DVec2) {
        let (p1, t1, d21) = curve_tool::d2(&self.curve1, x[0]);
        let (p2, t2, d22) = curve_tool::d2(&self.curve2, x[1]);
        (p1, p2, t1, t2, d21, d22)
    }
}

impl FunctionSetWithDerivatives for FunctionTanCuCu {
    fn nb_variables(&self) -> usize {
        2
    }
    fn nb_equations(&self) -> usize {
        2
    }

    /// OCCT Value(X, Fval).
    fn value(&mut self, x: &[f64], fval: &mut [f64]) -> bool {
        let (point1, point2, vect11, vect21, _, _) = self.init_derivative(x);
        let norme_d11 = vect11.length();
        let norme_d21 = vect21.length();
        let the_direction = point2 - point1;
        let squaredir = the_direction.dot(the_direction);
        fval[0] = cross2d(the_direction, vect11) / (norme_d11 * squaredir);
        fval[1] = cross2d(vect11, vect21) / (norme_d11 * norme_d21);
        true
    }

    /// OCCT Derivatives(X, Deriv).
    fn derivatives(&mut self, x: &[f64], deriv: &mut [Vec<f64>]) -> bool {
        let (point1, point2, vect11, vect21, vect12, vect22) = self.init_derivative(x);
        let norme_d11 = vect11.length();
        let norme_d21 = vect21.length();
        let the_direction = point2 - point1;
        let squaredir = the_direction.dot(the_direction);
        deriv[0][0] = cross2d(the_direction, vect12) / (norme_d11 * squaredir)
            + (cross2d(the_direction, vect11) * norme_d11 * norme_d11
                * vect11.dot(the_direction))
                / (norme_d11 * norme_d11 * norme_d11 * squaredir * squaredir * squaredir);
        deriv[0][1] = cross2d(vect21, vect11) / (norme_d11 * squaredir)
            - (cross2d(the_direction, vect11) * norme_d11 * norme_d11
                * vect21.dot(the_direction))
                / (norme_d11 * norme_d11 * norme_d11 * squaredir * squaredir * squaredir);
        deriv[1][0] = cross2d(vect12, vect21) / (norme_d11 * norme_d21)
            - cross2d(vect11, vect21) * vect12.dot(vect11) * norme_d21 * norme_d21
                / (norme_d11 * norme_d11 * norme_d11 * norme_d21 * norme_d21 * norme_d21);
        deriv[1][1] = cross2d(vect11, vect22) / (norme_d11 * norme_d21)
            - cross2d(vect11, vect21) * vect22.dot(vect21) * norme_d11 * norme_d11
                / (norme_d11 * norme_d11 * norme_d11 * norme_d21 * norme_d21 * norme_d21);
        true
    }

    /// OCCT Values(X, Fval, Deriv).
    fn values(&mut self, x: &[f64], fval: &mut [f64], deriv: &mut [Vec<f64>]) -> bool {
        self.value(x, fval);
        self.derivatives(x, deriv);
        true
    }
}

/// 2D cross product (OCCT gp_Vec2d::Crossed).
fn cross2d(a: DVec2, b: DVec2) -> f64 {
    a.x * b.y - a.y * b.x
}

/// Signed angle from a to b in [-π, π] (OCCT gp_Vec2d::Angle).
fn angle2d(a: DVec2, b: DVec2) -> f64 {
    (a.x * b.y - a.y * b.x).atan2(a.dot(b))
}

/// OCCT math_FunctionRoot (math_FunctionRoot.cxx L73-119) — the bounded
/// 1-variable root via the FunctionSetRoot (1 equation).
fn function_root_1d(
    f: &mut dyn FnMut(f64) -> f64,
    f_deriv: &mut dyn FnMut(f64) -> f64,
    guess: f64,
    tolerance: f64,
    a: f64,
    b: f64,
    nb_iterations: i32,
) -> Option<(f64, f64, f64)> {
    struct Adapter<'a> {
        f: &'a mut dyn FnMut(f64) -> f64,
        f_deriv: &'a mut dyn FnMut(f64) -> f64,
    }
    impl<'a> FunctionSetWithDerivatives for Adapter<'a> {
        fn nb_variables(&self) -> usize {
            1
        }
        fn nb_equations(&self) -> usize {
            1
        }
        fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
            f[0] = (self.f)(x[0]);
            true
        }
        fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
            df[0][0] = (self.f_deriv)(x[0]);
            true
        }
        fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
            f[0] = (self.f)(x[0]);
            df[0][0] = (self.f_deriv)(x[0]);
            true
        }
    }
    let mut adapter = Adapter {
        f,
        f_deriv,
    };
    let mut sol = FunctionSetRoot::new(&adapter, &[tolerance], nb_iterations);
    sol.perform(&mut adapter, &[guess], &[a], &[b], false);
    if !sol.is_done() {
        return None;
    }
    let root = sol.root()[0];
    let deriv = sol.derivative();
    let error = f(root);
    Some((root, deriv, error))
}

/// OCCT Geom2dGcc_Lin2d2TanIter (Geom2dGcc_Lin2d2TanIter.cxx, whole) — the
/// iterative tangent-line solver for (curve, point) and (curve, curve).
#[derive(Debug, Clone)]
pub struct Lin2d2TanIter {
    done: bool,
    line: (DVec2, DVec2),
    qualifier1: GccEntPosition,
    qualifier2: GccEntPosition,
    pnt_tg1_sol: DVec2,
    pnt_tg2_sol: DVec2,
    par1_sol: f64,
    par2_sol: f64,
    par_arg1: f64,
    par_arg2: f64,
}

impl Lin2d2TanIter {
    fn new_empty() -> Self {
        Lin2d2TanIter {
            done: false,
            line: (DVec2::ZERO, DVec2::X),
            qualifier1: GccEntPosition::NoQualifier,
            qualifier2: GccEntPosition::NoQualifier,
            pnt_tg1_sol: DVec2::ZERO,
            pnt_tg2_sol: DVec2::ZERO,
            par1_sol: 0.0,
            par2_sol: 0.0,
            par_arg1: 0.0,
            par_arg2: 0.0,
        }
    }

    /// OCCT Geom2dGcc_Lin2d2TanIter(Qualified1, ThePoint, Param1, Tolang)
    /// (L218-282) — tangent line through the point to the curve.
    pub fn curve_point(q1: &QCurve, the_point: DVec2, param1: f64, tolang: f64) -> Self {
        let mut r = Lin2d2TanIter::new_empty();
        if !(q1.is_enclosed() || q1.is_enclosing() || q1.is_outside() || q1.is_unqualified()) {
            panic!("GccEnt_BadQualifier");
        }
        let cu1 = q1.qualified().clone();
        let u1 = curve_tool::first_parameter(&cu1);
        let u2 = curve_tool::last_parameter(&cu1);

        let func = FunctionTanCuPnt::new(cu1.clone(), the_point);
        let mut fv = |x: f64| func.value(x);
        let mut fd = |x: f64| func.derivative(x);
        let eps_x = curve_tool::eps_x(&cu1, tolang.abs());
        if let Some((u_sol, _, norm)) = function_root_1d(&mut fv, &mut fd, param1, eps_x, u1, u2, 100) {
            if norm.abs() < tolang {
                let (origine, vect1, vect2) = curve_tool::d2(&cu1, u_sol);
                let vdir = the_point - origine;
                let mut sign1 = vect1.dot(vdir);
                let sign2 = cross2d(vect2, vdir);
                let accept = q1.is_unqualified()
                    || (q1.is_enclosing() && ((sign1 >= 0.0 && sign2 <= 0.0) || (sign1 <= 0.0 && sign2 <= 0.0)))
                    || (q1.is_outside() && sign1 <= 0.0 && sign2 >= 0.0)
                    || (q1.is_enclosed() && sign1 >= 0.0 && sign2 >= 0.0);
                if accept {
                    r.done = true;
                    r.line = (origine, vdir.normalize_or_zero());
                    r.qualifier1 = q1.qualifier();
                    r.qualifier2 = GccEntPosition::NoQualifier;
                    r.pnt_tg1_sol = origine;
                    r.pnt_tg2_sol = the_point;
                    r.par_arg1 = u_sol;
                    r.par1_sol = 0.0;
                    r.par_arg2 = the_point.distance(origine);
                    r.par2_sol = 0.0;
                }
            }
        }
        r
    }

    /// OCCT Geom2dGcc_Lin2d2TanIter(Qualified1, Qualified2, Param1, Param2,
    /// Tolang) (L137-216) — common tangent line to two curves.
    pub fn curve_curve(q1: &QCurve, q2: &QCurve, param1: f64, param2: f64, tolang: f64) -> Self {
        let mut r = Lin2d2TanIter::new_empty();
        if !(q1.is_enclosed() || q1.is_enclosing() || q1.is_outside() || q1.is_unqualified())
            || !(q2.is_enclosed() || q2.is_enclosing() || q2.is_outside() || q2.is_unqualified())
        {
            panic!("GccEnt_BadQualifier");
        }
        let cu1 = q1.qualified().clone();
        let cu2 = q2.qualified().clone();

        let mut func = FunctionTanCuCu::new(cu1.clone(), cu2.clone());
        let umin = [
            curve_tool::first_parameter(&cu1),
            curve_tool::first_parameter(&cu2),
        ];
        let umax = [
            curve_tool::last_parameter(&cu1),
            curve_tool::last_parameter(&cu2),
        ];
        let ufirst = [param1, param2];
        let tol = [
            curve_tool::eps_x(&cu1, tolang.abs()),
            curve_tool::eps_x(&cu2, tolang.abs()),
        ];

        let mut root = FunctionSetRoot::new(&func, &tol, 100);
        root.perform(&mut func, &ufirst, &umin, &umax, false);
        if root.is_done() {
            let ufirst = root.root();
            let mut norm = [0.0; 2];
            func.value(&ufirst, &mut norm);
            if norm[0].abs() < tolang && norm[1].abs() < tolang {
                // OCCT: D2 of both curves at the solution.
                let (point1, point2, vect11, vect21, vect12, vect22) =
                    curve_tool_d2_pair(&cu1, &cu2, ufirst[0], ufirst[1]);
                // OCCT Vec(point1 → point2).
                let vec = point2 - point1;
                let mut angle1 = angle2d(vec, vect12);
                let mut sign1 = vect11.dot(vec);
                let accept1 = q1.is_unqualified()
                    || (q1.is_enclosing() && angle1 >= 0.0)
                    || (q1.is_outside() && angle1 <= 0.0 && sign1 <= 0.0)
                    || (q1.is_enclosed() && angle1 <= 0.0 && sign1 >= 0.0);
                if accept1 {
                    angle1 = angle2d(vec, vect22);
                    sign1 = vect21.dot(vec);
                    let accept2 = q2.is_unqualified()
                        || (q2.is_enclosing() && angle1 >= 0.0)
                        || (q2.is_outside() && angle1 <= 0.0 && sign1 <= 0.0)
                        || (q2.is_enclosed() && angle1 <= 0.0 && sign1 >= 0.0);
                    if accept2 {
                        r.qualifier1 = q1.qualifier();
                        r.qualifier2 = q2.qualifier();
                        r.par_arg1 = ufirst[0];
                        r.par1_sol = 0.0;
                        r.pnt_tg1_sol = point1;
                        r.par_arg2 = ufirst[1];
                        r.pnt_tg2_sol = point2;
                        r.par2_sol = point2.distance(point1);
                        r.line = (point1, (point2 - point1).normalize_or_zero());
                        r.done = true;
                    }
                }
            }
        }
        r
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT ThisSolution() — the line (origin, direction).
    pub fn this_solution(&self) -> (DVec2, DVec2) {
        assert!(self.done, "StdFail_NotDone");
        self.line
    }

    /// OCCT WhichQualifier(Qualif1, Qualif2).
    pub fn which_qualifier(&self) -> (GccEntPosition, GccEntPosition) {
        assert!(self.done, "StdFail_NotDone");
        (self.qualifier1, self.qualifier2)
    }

    /// OCCT Tangency1(ParSol, ParArg, Pnt).
    pub fn tangency1(&self) -> (f64, f64, DVec2) {
        assert!(self.done, "StdFail_NotDone");
        (self.par1_sol, self.par_arg1, self.pnt_tg1_sol)
    }

    /// OCCT Tangency2(ParSol, ParArg, Pnt).
    pub fn tangency2(&self) -> (f64, f64, DVec2) {
        assert!(self.done, "StdFail_NotDone");
        (self.par2_sol, self.par_arg2, self.pnt_tg2_sol)
    }
}

/// D2 of two curves at two parameters.
fn curve_tool_d2_pair(
    c1: &Curve2d,
    c2: &Curve2d,
    u1: f64,
    u2: f64,
) -> (DVec2, DVec2, DVec2, DVec2, DVec2, DVec2) {
    let (p1, t1, d21) = curve_tool::d2(c1, u1);
    let (p2, t2, d22) = curve_tool::d2(c2, u2);
    (p1, p2, t1, t2, d21, d22)
}

/// OCCT Geom2dGcc_Lin2d2Tan (Geom2dGcc_Lin2d2Tan.cxx, whole) — tangent line to
/// a curve and a point / to two curves.
#[derive(Debug, Clone)]
pub struct Lin2d2Tan {
    well_done: bool,
    nbr_sol: usize,
    lines: Vec<(DVec2, DVec2)>,
    qualifiers1: Vec<GccEntPosition>,
    qualifiers2: Vec<GccEntPosition>,
    pnt_tg1: Vec<DVec2>,
    pnt_tg2: Vec<DVec2>,
    par1_sol: Vec<f64>,
    par2_sol: Vec<f64>,
    par_arg1: Vec<f64>,
    par_arg2: Vec<f64>,
}

impl Lin2d2Tan {
    fn new_with_capacity(cap: usize) -> Self {
        Lin2d2Tan {
            well_done: false,
            nbr_sol: 0,
            lines: Vec::with_capacity(cap),
            qualifiers1: Vec::with_capacity(cap),
            qualifiers2: Vec::with_capacity(cap),
            pnt_tg1: Vec::with_capacity(cap),
            pnt_tg2: Vec::with_capacity(cap),
            par1_sol: Vec::with_capacity(cap),
            par2_sol: Vec::with_capacity(cap),
            par_arg1: Vec::with_capacity(cap),
            par_arg2: Vec::with_capacity(cap),
        }
    }

    /// OCCT Geom2dGcc_Lin2d2Tan(Qualified1, ThePoint, Tolang) (L122-191) —
    /// the curve is sampled and the iterative solver runs from each start.
    pub fn curve_point(q1: &QualifiedCurve, the_point: DVec2, tolang: f64) -> Self {
        let mut r = Lin2d2Tan::new_with_capacity(2);
        let c1 = q1.qualified().clone();

        r.nbr_sol = 0;
        if let Curve2d::Circle(cir) = &c1 {
            // OCCT analytic case via GccAna_Lin2d2Tan (circle, point)
            // (GccAna_Lin2d2Tan.cxx L83-174).
            let circ = *cir;
            let r1 = circ.radius;
            let dist = the_point.distance(circ.center);
            if tolang.abs() < r1 - dist {
                r.well_done = true;
            } else if (dist - r1).abs() <= tolang.abs() {
                // The point lies on the circle: the tangent is perpendicular
                // to the radius (L118-129).
                let dir = (the_point - circ.center).normalize_or_zero();
                let line = (the_point, DVec2::new(-dir.y, dir.x));
                r.lines.push(line);
                r.qualifiers1.push(q1.qualifier());
                r.qualifiers2.push(GccEntPosition::NoQualifier);
                r.pnt_tg1.push(the_point);
                r.pnt_tg2.push(the_point);
                r.par_arg1.push(0.0);
                r.par_arg2.push(0.0);
                r.nbr_sol = 1;
                r.well_done = true;
            } else {
                // Two tangents from the point (L130-166).
                let d = dist - (dist * dist - r1 * r1).sqrt();
                let (mut signe, nbr_sol) = match q1.qualifier() {
                    GccEntPosition::Enclosing => (1.0, 1),
                    GccEntPosition::Outside => (-1.0, 1),
                    _ => (1.0, 2),
                };
                for _ in 0..nbr_sol {
                    // P1 = C1.Location().Rotated(ThePoint, asin(signe·R1/dist)).
                    let rel = circ.center - the_point;
                    let alpha = (signe * r1 / dist).asin();
                    let rotated = the_point
                        + DVec2::new(
                            rel.x * alpha.cos() - rel.y * alpha.sin(),
                            rel.x * alpha.sin() + rel.y * alpha.cos(),
                        );
                    let mut dir = (the_point - rotated).normalize_or_zero();
                    let p1 = rotated + d * dir;
                    dir = (the_point - p1).normalize_or_zero();
                    r.lines.push((p1, dir));
                    r.qualifiers1.push(q1.qualifier());
                    r.qualifiers2.push(GccEntPosition::NoQualifier);
                    r.pnt_tg1.push(p1);
                    r.pnt_tg2.push(the_point);
                    r.par_arg1.push(0.0);
                    r.par_arg2.push(0.0);
                    r.nbr_sol += 1;
                    signe = -signe;
                }
                r.well_done = true;
            }
        } else {
            let qc1 = QCurve::new(c1.clone(), q1.qualifier());
            let a_first_par = curve_tool::first_parameter(&c1);
            let a_last_par = curve_tool::last_parameter(&c1);
            let a_nb_samples = curve_tool::nb_samples(&c1);
            let a_step = (a_last_par - a_first_par) / a_nb_samples as f64;
            let mut param1 = a_first_par;

            let mut i = 0;
            while i <= a_nb_samples && r.nbr_sol < 2 {
                let lin = Lin2d2TanIter::curve_point(&qc1, the_point, param1, tolang);
                if lin.is_done() && r.add(&lin, tolang, &c1, None) {
                    r.nbr_sol += 1;
                }
                param1 += a_step;
                i += 1;
            }
            r.well_done = r.nbr_sol > 0;
        }
        r
    }

    /// OCCT Geom2dGcc_Lin2d2Tan(Qualified1, Qualified2, Tolang) (L32-120) —
    /// analytic circle-circle via GccAna_Lin2d2Tan, general curves via the
    /// 20×20 sampling + Lin2d2TanIter + Add.
    pub fn curve_curve(q1: &QualifiedCurve, q2: &QualifiedCurve, tolang: f64) -> Self {
        let mut r = Lin2d2Tan::new_with_capacity(4);
        let c1 = q1.qualified().clone();
        let c2 = q2.qualified().clone();

        r.nbr_sol = 0;
        if let (Curve2d::Circle(cir1), Curve2d::Circle(cir2)) = (&c1, &c2) {
            // OCCT GccAna_Lin2d2Tan (circle, circle) — common tangents.
            let (c1c, r1) = (cir1.center, cir1.radius);
            let (c2c, r2) = (cir2.center, cir2.radius);
            let dist = c2c.distance(c1c);
            if dist >= tolang.abs() {
                let mut lines: Vec<((DVec2, DVec2), DVec2, DVec2)> = Vec::new();
                // External tangents: angle = asin((r2-r1)/dist).
                let base = (c2c - c1c).normalize_or_zero();
                let perp = DVec2::new(-base.y, base.x);
                let add_ext = |lines: &mut Vec<((DVec2, DVec2), DVec2, DVec2)>,
                               sign: f64| {
                    let alpha = ((r2 - r1) / dist).asin();
                    let dir = base * alpha.cos() + perp * (sign * alpha.sin());
                    let n1 = DVec2::new(-dir.y, dir.x);
                    let p1 = c1c + n1 * r1;
                    let p2 = c2c + n1 * r2;
                    lines.push(((p1, dir), p1, p2));
                };
                add_ext(&mut lines, 1.0);
                add_ext(&mut lines, -1.0);
                // Internal tangents when the circles are separate.
                if r1 + r2 < dist {
                    let add_int = |lines: &mut Vec<((DVec2, DVec2), DVec2, DVec2)>,
                                   sign: f64| {
                        let alpha = ((r2 + r1) / dist).asin();
                        let dir = base * alpha.cos() + perp * (sign * alpha.sin());
                        let n1 = DVec2::new(-dir.y, dir.x);
                        let p1 = c1c + n1 * r1;
                        let p2 = c2c - n1 * r2;
                        lines.push(((p1, dir), p1, p2));
                    };
                    add_int(&mut lines, 1.0);
                    add_int(&mut lines, -1.0);
                }
                for (line, p1, p2) in lines {
                    r.lines.push(line);
                    r.qualifiers1.push(q1.qualifier());
                    r.qualifiers2.push(q2.qualifier());
                    r.pnt_tg1.push(p1);
                    r.pnt_tg2.push(p2);
                    r.par_arg1.push(0.0);
                    r.par_arg2.push(0.0);
                    r.nbr_sol += 1;
                }
                r.well_done = r.nbr_sol > 0;
            }
        } else {
            let qc1 = QCurve::new(c1.clone(), q1.qualifier());
            let qc2 = QCurve::new(c2.clone(), q2.qualifier());
            let a1_f_par = curve_tool::first_parameter(&c1);
            let a1_l_par = curve_tool::last_parameter(&c1);
            let a_nb_samples1 = curve_tool::nb_samples(&c1);
            let a_step1 = (a1_l_par - a1_f_par) / a_nb_samples1 as f64;
            let a2_f_par = curve_tool::first_parameter(&c2);
            let a2_l_par = curve_tool::last_parameter(&c2);
            let a_nb_samples2 = curve_tool::nb_samples(&c2);
            let a_step2 = (a2_l_par - a2_f_par) / a_nb_samples2 as f64;

            let mut param1 = a1_f_par;
            let mut i = 0;
            while i <= a_nb_samples1 && r.nbr_sol < 4 {
                let mut param2 = a2_f_par;
                let mut j = 0;
                while j <= a_nb_samples2 && r.nbr_sol < 4 {
                    let lin = Lin2d2TanIter::curve_curve(&qc1, &qc2, param1, param2, tolang);
                    if lin.is_done() && r.add(&lin, tolang, &c1, Some(&c2)) {
                        r.nbr_sol += 1;
                    }
                    param2 += a_step2;
                    j += 1;
                }
                param1 += a_step1;
                i += 1;
            }
            r.well_done = r.nbr_sol > 0;
        }
        r
    }

    /// OCCT Geom2dGcc_Lin2d2Tan::Add (L407-465) — dedup + tangency direction
    /// verification, then store the solution.
    fn add(
        &mut self,
        lin: &Lin2d2TanIter,
        the_tol: f64,
        c1: &Curve2d,
        c2: Option<&Curve2d>,
    ) -> bool {
        let (a_par1sol, a_par1arg, a_pnt1_sol) = lin.tangency1();
        let (a_par2sol, a_par2arg, a_pnt2_sol) = lin.tangency2();
        let (a_lin, a_lin_dir) = lin.this_solution();

        for i in 0..self.nbr_sol {
            if (a_par1arg - self.par_arg1[i]).abs() <= the_tol
                && (a_par2arg - self.par_arg2[i]).abs() <= the_tol
            {
                return false;
            }
        }

        let (a_point, a_vtan) = curve_tool::d1(c1, a_par1arg);
        let _ = a_point;
        if cross2d(a_lin_dir, a_vtan.normalize_or_zero()).abs() > the_tol {
            return false;
        }

        if let Some(c2) = c2 {
            let (a_point, a_vtan) = curve_tool::d1(c2, a_par2arg);
            let _ = a_point;
            if cross2d(a_lin_dir, a_vtan.normalize_or_zero()).abs() > the_tol {
                return false;
            }
        }

        let (q1, q2) = lin.which_qualifier();
        self.lines.push((a_lin, a_lin_dir));
        self.qualifiers1.push(q1);
        self.qualifiers2.push(q2);
        self.pnt_tg1.push(a_pnt1_sol);
        self.pnt_tg2.push(a_pnt2_sol);
        self.par1_sol.push(a_par1sol);
        self.par2_sol.push(a_par2sol);
        self.par_arg1.push(a_par1arg);
        self.par_arg2.push(a_par2arg);
        true
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.well_done
    }

    /// OCCT NbSolutions().
    pub fn nb_solutions(&self) -> usize {
        self.nbr_sol
    }
}
