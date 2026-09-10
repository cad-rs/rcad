//! OCCT GeomFill_FunctionGuide (TKGeomAlgo/GeomFill) — 1:1 port of
//! GeomFill_FunctionGuide.hxx (members) + GeomFill_FunctionGuide.cxx (whole
//! file L45-273).  The Deriv2T / DerivTX / Deriv2X bodies are commented out
//! in the OCCT source (L280-363) and are not compiled there either.
//!
//! Architecture differences:
//! - The class derives `math_FunctionSetWithDerivatives`; the rcad form
//!   implements [`FunctionSetWithDerivatives`] on the same struct
//!   (3 variables (w, u, v) / 3 equations).
//! - `TheGuide` / `TheLaw` are OCCT handles; here they are shared borrows
//!   (the function never owns or mutates them).
//! - `Geom_SurfaceOfRevolution` maps to `Surface3::Revolution` (u = the
//!   rotation angle, v = the basis-curve parameter — the same parameter
//!   convention as the OCCT surface).
//! - `gp_Trsf::SetTransformation(Ax3, Ax3)` maps to the local
//!   [`gp_trsf_set_transformation_between`] re-host (pure math; the
//!   brep_fill_evolved_d precedent).
//! - `TheCurve->Transform(Transfo)` maps to [`curve_transformed_by_trsf`]
//!   (BSpline poles / trimmed basis; the kernel transform_curve mapping).

use glam::{DAffine3, DMat3, DVec3};

use rcad_kernel::geom::{
    transform_curve, BSplineCurve3, Curve3, CurveEval, RevolutionSurface, Surface3, SurfaceEval,
    TrimmedCurve3,
};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;
use rcad_kernel::math::gp::{Ax1, Ax3};

use super::section_law::SectionLaw;
use super::trihedron_law::{curve_first_parameter, curve_last_parameter};

// ---------------------------------------------------------------------------
// Pure-math gp re-hosts (shared with location_guide.rs)
// ---------------------------------------------------------------------------

/// The gp_Ax3 frame affine (local -> world).
fn ax3_frame(a: &Ax3) -> DAffine3 {
    let m = DMat3::from_cols(a.x_direction, a.y_direction, a.axis.direction);
    let mut f = DAffine3::from_mat3(m);
    f.translation = a.axis.location;
    f
}

/// OCCT gp_Trsf::SetTransformation(FromT1, ToT2) — the transformation from
/// the T1 coordinate system to the T2 one (p' = Frame(T2)^-1 * Frame(T1) * p)
/// (pure-math re-host).
pub(crate) fn gp_trsf_set_transformation_between(from_t1: &Ax3, to_t2: &Ax3) -> DAffine3 {
    ax3_frame(to_t2).inverse() * ax3_frame(from_t1)
}

/// OCCT Geom_Curve::Transform(Trsf) over the rcad Curve3 forms consumed by
/// the guide functions — BSpline poles and the Trimmed basis curve
/// (Geom_TrimmedCurve::Transform keeps the trim range and transforms the
/// basis curve).
pub(crate) fn curve_transformed_by_trsf(c: &Curve3, trsf: &DAffine3) -> Curve3 {
    match c {
        Curve3::Trimmed(tc) => Curve3::Trimmed(TrimmedCurve3::new(
            curve_transformed_by_trsf(&tc.curve, trsf),
            tc.first,
            tc.last,
        )),
        other => transform_curve(other, trsf),
    }
}

/// OCCT GeomFill_FunctionGuide (GeomFill_FunctionGuide.hxx L86-113).
pub struct FunctionGuide<'a> {
    /// OCCT handle(Adaptor3d_Curve) TheGuide.
    the_guide: &'a Curve3,
    /// OCCT handle(GeomFill_SectionLaw) TheLaw.
    the_law: &'a dyn SectionLaw,
    /// OCCT bool isconst.
    isconst: bool,
    /// OCCT handle(Geom_Curve) TheCurve — the generatrix (set by SetParam).
    the_curve: Option<Curve3>,
    /// OCCT handle(Geom_Curve) TheConst — the constant section.
    the_const: Option<Curve3>,
    /// OCCT handle(Geom_Surface) TheSurface — the revolution surface.
    the_surface: Option<Surface3>,
    /// OCCT double First.
    first: f64,
    /// OCCT double Last.
    last: f64,
    /// OCCT double TheUonS.
    the_uon_s: f64,
    /// OCCT gp_XYZ Centre.
    centre: DVec3,
    /// OCCT gp_XYZ Dir.
    dir: DVec3,
}

impl<'a> FunctionGuide<'a> {
    /// OCCT GeomFill_FunctionGuide::GeomFill_FunctionGuide (L45-69).
    pub fn new(s: &'a dyn SectionLaw, guide: &'a Curve3, param_on_law: f64) -> Self {
        let mut function = FunctionGuide {
            the_guide: guide,
            the_law: s,
            isconst: false,
            the_curve: None,
            the_const: None,
            the_surface: None,
            first: 0.0,
            last: 0.0,
            the_uon_s: param_on_law,
            centre: DVec3::ZERO,
            dir: DVec3::ZERO,
        };
        // OCCT: double Tol = Precision::Confusion(); IsConstant(Tol) takes
        // the Error reference and may overwrite it.
        let mut tol = 1.0e-7;
        if function.the_law.is_constant(&mut tol) {
            // OCCT: isconst = true; TheConst = TheLaw->ConstantSection();
            // First = TheConst->FirstParameter(); Last = ...
            function.isconst = true;
            let the_const = function.the_law.constant_section();
            function.first = curve_first_parameter(&the_const);
            function.last = curve_last_parameter(&the_const);
            function.the_const = Some(the_const);
        } else {
            // OCCT: isconst = false; TheConst.Nullify().
            function.isconst = false;
            function.the_const = None;
        }
        // OCCT: TheCurve.Nullify().
        function.the_curve = None;
        function
    }

    /// OCCT GeomFill_FunctionGuide::SetParam (L76-121) — initialisation of
    /// the revolution surface (the Param argument is unused in the OCCT
    /// body).
    pub fn set_param(&mut self, _param: f64, c: DVec3, d: DVec3, dx: DVec3) {
        // OCCT: Centre = C.XYZ(); Dir = D.
        self.centre = c;
        self.dir = d;

        // repere fixe: gp_Ax3 Rep(gp::Origin(), gp::DZ(), gp::DX()).
        let rep = Ax3::new();

        // calculer transfo entre triedre et Oxyz:
        // gp_Dir B2 = DX; gp_Ax3 RepTriedre(C, D, B2);
        // gp_Trsf Transfo; Transfo.SetTransformation(RepTriedre, Rep).
        let b2 = dx.normalize_or_zero();
        let rep_triedre = Ax3::from_pnt_n_vx(c, d.normalize_or_zero(), b2);
        let transfo = gp_trsf_set_transformation_between(&rep_triedre, &rep);

        if self.isconst {
            // OCCT: TheCurve = new Geom_TrimmedCurve(TheConst->Copy(),
            // First, Last).
            let the_const = self.the_const.as_ref().expect("null TheConst");
            self.the_curve = Some(Curve3::Trimmed(TrimmedCurve3::new(
                the_const.clone(),
                self.first,
                self.last,
            )));
        } else {
            // OCCT: TheLaw->SectionShape(NbPoles, NbKnots, Deg);
            // Mults/Knots/D0(TheUonS, Poles, Weights); new
            // Geom_BSplineCurve(Poles, [Weights,] Knots, Mult, Deg,
            // TheLaw->IsUPeriodic()).
            let mut nb_poles = 0usize;
            let mut nb_knots = 0usize;
            let mut deg = 0usize;
            self.the_law.section_shape(&mut nb_poles, &mut nb_knots, &mut deg);
            let mut mults = vec![0i32; nb_knots];
            self.the_law.mults(&mut mults);
            let mut knots = vec![0.0f64; nb_knots];
            self.the_law.knots(&mut knots);
            let mut poles = vec![DVec3::ZERO; nb_poles];
            let mut weights = vec![1.0f64; nb_poles];
            self.the_law.d0(self.the_uon_s, &mut poles, &mut weights);
            let mut bs = BSplineCurve3::from_knots_mults(deg, knots, mults, poles);
            if self.the_law.is_rational() {
                bs.weights = weights;
            }
            bs.is_periodic = self.the_law.is_u_periodic();
            self.the_curve = Some(Curve3::BSpline(bs));
        }

        // OCCT: gp_Ax1 Axe(C, Dir); TheCurve->Transform(Transfo);
        // TheSurface = new Geom_SurfaceOfRevolution(TheCurve, Axe).
        let the_curve = self.the_curve.as_ref().expect("null TheCurve");
        let transformed = curve_transformed_by_trsf(the_curve, &transfo);
        let axe = Ax1::new(c, self.dir);
        self.the_surface = Some(Surface3::Revolution(RevolutionSurface {
            profile: Box::new(transformed),
            axis_origin: axe.location,
            axis_dir: axe.direction,
        }));
    }

    /// OCCT GeomFill_FunctionGuide::Value (L143-155) — G(w) - S(u, v).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        // TheGuide->D0(X(1), P).
        let p = self.the_guide.point_at(x[0]);
        // TheSurface->D0(X(2), X(3), P1).
        let the_surface = self.the_surface.as_ref().expect("null TheSurface");
        let p1 = the_surface.point_at(x[1], x[2]);

        f[0] = p.x - p1.x;
        f[1] = p.y - p1.y;
        f[2] = p.z - p1.z;

        true
    }

    /// OCCT GeomFill_FunctionGuide::Derivatives (L161-178).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        // TheGuide->D1(X(1), P, DP).
        let p = self.the_guide.point_at(x[0]);
        let dp = self.the_guide.derivative_at(x[0]);
        // TheSurface->D1(X(2), X(3), P1, DP1U, DP1V).
        let the_surface = self.the_surface.as_ref().expect("null TheSurface");
        let (_p1, dp1u, dp1v) = the_surface.derivatives(x[1], x[2]);

        for i in 0..3 {
            d[i][0] = dp[i];
            d[i][1] = -dp1u[i];
            d[i][2] = -dp1v[i];
        }

        true
    }

    /// OCCT GeomFill_FunctionGuide::Values (L184-203).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        // TheGuide->D1(X(1), P, DP) — derivee de la generatrice.
        let p = self.the_guide.point_at(x[0]);
        let dp = self.the_guide.derivative_at(x[0]);
        // TheSurface->D1(X(2), X(3), P1, DP1U, DP1V) — derivee de la new
        // surface.
        let the_surface = self.the_surface.as_ref().expect("null TheSurface");
        let (p1, dp1u, dp1v) = the_surface.derivatives(x[1], x[2]);

        for i in 0..3 {
            f[i] = p[i] - p1[i];
            d[i][0] = dp[i];
            d[i][1] = -dp1u[i];
            d[i][2] = -dp1v[i];
        }

        true
    }

    /// OCCT GeomFill_FunctionGuide::DerivT (L209-225) — the first derivative
    /// from t.
    pub fn deriv_t(&mut self, x: &[f64], d_centre: DVec3, d_dir: DVec3, f: &mut [f64]) -> bool {
        // DSDT(X(2), X(3), DCentre, DDir, DS).
        let mut ds = DVec3::ZERO;
        self.dsdt(x[1], x[2], d_centre, d_dir, &mut ds);

        // TheCurve->D0(X(1), P).
        let the_curve = self.the_curve.as_ref().expect("null TheCurve");
        let p = the_curve.point_at(x[0]);

        f[0] = p.x - ds.x;
        f[1] = p.y - ds.y;
        f[2] = p.z - ds.z;

        true
    }

    /// OCCT GeomFill_FunctionGuide::DSDT (L231-273) — the derivative of the
    /// revolution surface with respect to t, at (U, V).
    pub fn dsdt(&self, u: f64, v: f64, dc: DVec3, ddir: DVec3, ds: &mut DVec3) {
        // TheCurve->D0(V, Pc) — Q(v).
        let the_curve = self.the_curve.as_ref().expect("null TheCurve");
        let pc = the_curve.point_at(v);
        let mut q = pc;
        let mut dq = DVec3::ZERO;
        if !self.isconst {
            // OCCT prints "Not implemented" for the moving (non-constant)
            // section and keeps DQ = 0.
        }

        // Q.Subtract(Centre) — CQ; DQ -= DC.
        q -= self.centre;
        dq -= dc;

        // DVcrossCQ.SetLinearForm(DDir.Crossed(Q), Dir.Crossed(DQ)) —
        // Vdir^CQ; DVcrossCQ.Multiply(sin(U)).
        let mut dvcrosscq = ddir.cross(q) + self.dir.cross(dq);
        dvcrosscq *= u.sin();

        // DVdotCQ.SetLinearForm(DDir.Dot(Q) + Dir.Dot(DQ), Dir,
        // Dir.Dot(Q), DDir) — (CQ.Vdir)(1-cos(U))Vdir.
        let mut dvdotcq = (ddir.dot(q) + self.dir.dot(dq)) * self.dir + self.dir.dot(q) * ddir;
        // DVdotCQ.Add(DVcrossCQ) — addition des composantes.
        dvdotcq += dvcrosscq;

        let cos_u = u.cos();
        dq *= cos_u;
        dq += dvdotcq;
        dq += dc;
        *ds = dq;
    }
}

impl FunctionSetWithDerivatives for FunctionGuide<'_> {
    /// OCCT NbVariables (L127-130) — (w, u, v).
    fn nb_variables(&self) -> usize {
        3
    }

    /// OCCT NbEquations (L134-137).
    fn nb_equations(&self) -> usize {
        3
    }

    /// OCCT Value.
    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        FunctionGuide::value(self, x, f)
    }

    /// OCCT Derivatives.
    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        FunctionGuide::derivatives(self, x, df)
    }

    /// OCCT Values.
    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        FunctionGuide::values(self, x, f, df)
    }
}
