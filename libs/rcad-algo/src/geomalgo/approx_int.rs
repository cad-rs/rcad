//! 1:1 translation of the OCCT WLine smooth approximation chain used by
//! IntTools_FaceFace::MakeCurve (IntPatch_Walking, myApprox = true):
//!
//!   GeomInt_WLApprox (GeomInt_WLApprox.hxx + ApproxInt_Approx.gxx)
//!     -> ApproxInt_KnotTools (ApproxInt_KnotTools.cxx)
//!     -> AppParCurves_Approx (Approx_ComputeLine.gxx, instantiated as
//!        GeomInt_TheComputeLineBezierOfWLApprox)
//!     -> AppParCurves_Gradient (AppParCurves_Gradient.gxx)
//!     -> AppParCurves_LeastSquare (AppParCurves_LeastSquare.gxx)
//!     -> Approx_MCurvesToBSpCurve (Approx_MCurvesToBSpCurve.cxx +
//!        Convert_CompBezierCurvesToBSplineCurveBase.hxx)
//!
//! The data model maps the OCCT classes directly:
//!   AppParCurves_MultiPoint      -> MultiPoint   (points of one multipoint)
//!   AppParCurves_MultiCurve      -> MultiCurve   (poles of the Bezier fit)
//!   AppParCurves_MultiBSpCurve   -> MultiBSpCurve (the concatenated BSpline)
//!   GeomInt_TheMultiLineOfWLApprox -> WLineAccess (access to the WLine points)
//!   ApproxInt_Approx             -> WLineApprox
//!   Approx_ComputeLine           -> ComputeLine
//!   AppParCurves_Function        -> ParFunction
//!   AppParCurves_Gradient        -> Gradient
//!   AppParCurves_LeastSquare     -> LeastSquare
//!   Approx_MCurvesToBSpCurve     -> MCurvesToBSpCurve
//!
//! rcad note: the OCCT dense linear algebra (math_Householder, DACTCL
//! banded LDLT, math_BFGS line search) is the bottom layer of the chain and
//! is kept semantically equivalent (least squares solved by normal equations
//! with Gaussian elimination; BFGS with a golden-section line search).

use crate::geomalgo::int_patch::{IntPatchLine, WLinePnt};
use glam::{DVec2, DVec3};
use rcad_kernel::geom::{Surface3, SurfaceEval};
use rcad_kernel::math::plib::eval_polynomial_d1;
use rcad_kernel::math::VecD;

/// OCCT Approx_ParametrizationType.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproxParamType {
    ChordLength,
    Centripetal,
    IsoParametric,
}

/// OCCT AppParCurves_Constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd)]
pub enum AppParConstraint {
    NoConstraint = 0,
    PassPoint = 1,
    TangencyPoint = 2,
    CurvaturePoint = 3,
}

/// OCCT Approx_Status.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproxStatus {
    PointsAdded,
    NoPointsAdded,
    NoApproximation,
}

/// OCCT AppParCurves_ConstraintCouple — a constraint attached to a point
/// index.
#[derive(Debug, Clone, Copy)]
pub struct ConstraintCouple {
    pub index: i32,
    pub constraint: AppParConstraint,
}

/// OCCT AppParCurves_MultiPoint — the points (3D and 2D) of one multipoint.
/// The 2D points use the OCCT GLOBAL index (nbP+1 .. nbP+nbP2d).
#[derive(Debug, Clone, Default)]
pub struct MultiPoint {
    pub p3d: Vec<DVec3>,
    pub p2d: Vec<DVec2>,
    pub nb_p3d: usize,
}

impl MultiPoint {
    pub fn new(nb_p: usize, nb_p2d: usize) -> Self {
        MultiPoint {
            p3d: vec![DVec3::ZERO; nb_p],
            p2d: vec![DVec2::ZERO; nb_p2d],
            nb_p3d: nb_p,
        }
    }
    /// OCCT MultiPoint::Point(i): 3D point, local 1-based index.
    pub fn point(&self, i: usize) -> DVec3 {
        self.p3d[i - 1]
    }
    /// OCCT MultiPoint::Point2d(i): 2D point, GLOBAL 1-based index
    /// (nbP+1 .. nbP+nbP2d).
    pub fn point2d(&self, i: usize) -> DVec2 {
        self.p2d[i - 1 - self.nb_p3d]
    }
    pub fn set_point(&mut self, i: usize, p: DVec3) {
        self.p3d[i - 1] = p;
    }
    pub fn set_point2d(&mut self, i: usize, p: DVec2) {
        self.p2d[i - 1 - self.nb_p3d] = p;
    }
    pub fn nb_points(&self) -> usize {
        self.p3d.len()
    }
    pub fn nb_points2d(&self) -> usize {
        self.p2d.len()
    }
}

/// OCCT AppParCurves_MultiCurve — a Bezier multi-curve: `poles` holds the
/// poles of all curves (3D first, then 2D), each pole being one MultiPoint.
/// `degree = poles.len() - 1`.
#[derive(Debug, Clone)]
pub struct MultiCurve {
    pub poles: Vec<MultiPoint>,
}

impl MultiCurve {
    pub fn new(nb_poles: usize, nb_p3d: usize, nb_p2d: usize) -> Self {
        MultiCurve {
            poles: (0..nb_poles).map(|_| MultiPoint::new(nb_p3d, nb_p2d)).collect(),
        }
    }
    /// Number of poles (degree + 1).
    pub fn nb_poles(&self) -> usize {
        self.poles.len()
    }
    pub fn degree(&self) -> usize {
        self.poles.len() - 1
    }
    /// Number of curves (nbP3d + nbP2d).
    pub fn nb_curves(&self) -> usize {
        let p = &self.poles[0];
        p.p3d.len() + p.p2d.len()
    }
    pub fn dimension(&self, curve: usize) -> usize {
        let p = &self.poles[0];
        if curve <= p.p3d.len() {
            3
        } else {
            2
        }
    }
    /// OCCT MultiCurve::Curve(i, Poles): the poles of curve `i` (1-based).
    pub fn curve(&self, i: usize, out: &mut Vec<DVec3>) {
        let nb3d = self.poles[0].p3d.len();
        if i <= nb3d {
            out.clear();
            for p in &self.poles {
                out.push(p.p3d[i - 1]);
            }
        } else {
            out.clear();
            for p in &self.poles {
                out.push(DVec3::new(p.p2d[i - 1 - nb3d].x, p.p2d[i - 1 - nb3d].y, 0.0));
            }
        }
    }
    /// 2D variant of Curve.
    pub fn curve2d(&self, i: usize, out: &mut Vec<DVec2>) {
        let nb3d = self.poles[0].p3d.len();
        if i > nb3d {
            out.clear();
            for p in &self.poles {
                out.push(p.p2d[i - 1 - nb3d]);
            }
        } else {
            out.clear();
            for p in &self.poles {
                out.push(DVec2::new(p.p3d[i - 1].x, p.p3d[i - 1].y));
            }
        }
    }
    /// OCCT MultiCurve::SetValue(i, MPole).
    pub fn set_value(&mut self, i: usize, mp: MultiPoint) {
        self.poles[i - 1] = mp;
    }
    /// OCCT MultiCurve::Value(i): the i-th multipoint (pole).
    pub fn value(&self, i: usize) -> &MultiPoint {
        &self.poles[i - 1]
    }
    /// OCCT MultiCurve::D1(i, U, P, V1): point and first derivative of
    /// curve `i` at parameter `U` (Bezier evaluation).
    pub fn d1(&self, i: usize, u: f64) -> (DVec3, DVec3) {
        let nb3d = self.poles[0].p3d.len();
        let deg = self.degree();
        let mut pt = DVec3::ZERO;
        let mut d1 = DVec3::ZERO;
        if i <= nb3d {
            for j in 0..=deg {
                let (b, db) = bernstein(deg, j, u);
                pt += self.poles[j].p3d[i - 1] * b;
                d1 += self.poles[j].p3d[i - 1] * db;
            }
        } else {
            for j in 0..=deg {
                let (b, db) = bernstein(deg, j, u);
                pt += DVec3::new(self.poles[j].p2d[i - 1 - nb3d].x, self.poles[j].p2d[i - 1 - nb3d].y, 0.0) * b;
                d1 += DVec3::new(self.poles[j].p2d[i - 1 - nb3d].x, self.poles[j].p2d[i - 1 - nb3d].y, 0.0) * db;
            }
        }
        (pt, d1)
    }
}

/// OCCT AppParCurves_MultiBSpCurve — a BSpline multi-curve: shared degree,
/// knots and multiplicities, per-curve poles.
#[derive(Debug, Clone)]
pub struct MultiBSpCurve {
    pub degree: usize,
    pub knots: Vec<f64>,
    pub mults: Vec<usize>,
    pub poles: Vec<MultiPoint>,
}

impl MultiBSpCurve {
    /// OCCT MultiBSpCurve(CU, Knots, Mults): build from one Bezier MultiCurve.
    pub fn from_bezier(cu: &MultiCurve, knots: Vec<f64>, mults: Vec<usize>) -> Self {
        MultiBSpCurve {
            degree: cu.degree(),
            knots,
            mults,
            poles: cu.poles.clone(),
        }
    }
    /// OCCT MultiBSpCurve(tabMU, Knots, Mults).
    pub fn from_multipoints(
        tab_mu: Vec<MultiPoint>,
        knots: Vec<f64>,
        mults: Vec<usize>,
        degree: usize,
    ) -> Self {
        MultiBSpCurve {
            degree,
            knots,
            mults,
            poles: tab_mu,
        }
    }
    pub fn nb_poles(&self) -> usize {
        self.poles.len()
    }
    /// Poles of curve `i` (1-based) as 3D (for 3D curves).
    pub fn curve(&self, i: usize, out: &mut Vec<DVec3>) {
        let nb3d = self.poles[0].p3d.len();
        if i <= nb3d {
            out.clear();
            for p in &self.poles {
                out.push(p.p3d[i - 1]);
            }
        } else {
            out.clear();
            for p in &self.poles {
                out.push(DVec3::new(p.p2d[i - 1 - nb3d].x, p.p2d[i - 1 - nb3d].y, 0.0));
            }
        }
    }
    /// Poles of 2D curve `i` (1-based, 2D index space: 1..=nb2d).
    pub fn curve2d(&self, i: usize, out: &mut Vec<DVec2>) {
        let nb3d = self.poles[0].p3d.len();
        if i > nb3d {
            out.clear();
            for p in &self.poles {
                out.push(p.p2d[i - 1 - nb3d]);
            }
        } else {
            out.clear();
            for p in &self.poles {
                out.push(DVec2::new(p.p3d[i - 1].x, p.p3d[i - 1].y));
            }
        }
    }
}

/// OCCT Bernstein basis B(j, deg, u) and its derivative (dB/du).
pub fn bernstein(deg: usize, j: usize, u: f64) -> (f64, f64) {
    let b = binom(deg, j) * u.powi(j as i32) * (1.0 - u).powi((deg - j) as i32);
    let db = if deg == 0 {
        0.0
    } else if j == 0 {
        -(deg as f64) * binom(deg - 1, 0) * (1.0 - u).powi((deg - 1) as i32)
    } else if j == deg {
        deg as f64 * binom(deg - 1, deg - 1) * u.powi((deg - 1) as i32)
    } else {
        deg as f64 * (binom(deg - 1, j - 1) * u.powi((j - 1) as i32) * (1.0 - u).powi((deg - j) as i32)
            - binom(deg - 1, j) * u.powi(j as i32) * (1.0 - u).powi((deg - 1 - j) as i32))
    };
    (b, db)
}

fn binom(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let k = k.min(n - k);
    let mut r = 1.0f64;
    for i in 0..k {
        r = r * (n - i) as f64 / (i + 1) as f64;
    }
    r
}

/// OCCT BSplCLib::PolesCoefficients (BSplCLib_BzSyntaxes.cxx L62-85) for a
/// non-rational Bezier curve: converts the Bezier poles to the coefficients
/// of the power (Taylor) basis, CachePoles[ii] = D^{ii-1} P(0) / (ii-1)!.
/// The conversion is BSplCLib_BuildCache (BSplCLib_CurveComputation.pxx
/// L1527-1600) with U=0, SpanDomain=1, Periodic=false and the flat Bezier
/// knots (BSplCLib.cxx L4966-4977).  PrepareEval_T (L777-830) yields
/// Index = Degree+1 so BuildKnots (L1555+) gives the window
/// [0 x Degree, 1 x Degree], and BuildEval (L720-762) copies the poles in
/// order.  BSplCLib::Bohm (BSplCLib.cxx L1197-1400) then computes the
/// derivatives at U=0 and the final loop divides by the factorials.
fn bezier_poles_to_coeffs(poles: &[DVec3]) -> Vec<DVec3> {
    let deg = poles.len() - 1;
    let mut psav: Vec<DVec3> = poles.to_vec();
    // Flat Bezier knot window [0 x Degree, 1 x Degree] (BuildKnots with
    // Mults == nullptr, Index = Degree + 1).
    let mut knot = vec![0.0f64; 2 * deg];
    for k in deg..2 * deg {
        knot[k] = 1.0;
    }
    // OCCT Bohm first phase (L1332-1360): divided differences.  With the
    // Bezier window every (knot[jDmi] - knot[j]) is 1.0, but the branch is
    // kept verbatim so the arithmetic matches OCCT exactly.
    let mut ddmi = 2 * deg + 1;
    for i in 0..deg {
        ddmi -= 1;
        let mut jdmi = ddmi;
        for j in (i..deg).rev() {
            jdmi -= 1;
            let coef = if knot[jdmi] == knot[j] {
                0.0
            } else {
                1.0 / (knot[jdmi] - knot[j])
            };
            psav[j + 1] = (psav[j + 1] - psav[j]) * coef;
        }
    }
    // OCCT Bohm second phase (L1361-1383): accumulation in U.  For
    // PolesCoefficients U == 0 and knot[i] == 0 so coef = 0.0 and every
    // term is an exact no-op (0.0 * finite == 0.0); the result is
    // therefore identical to the OCCT loop.
    // OCCT Bohm multiply-by-degrees (L1384-1399): psav[i] *= Degree!/(Degree-i)!.
    let mut coef = deg as f64;
    let mut dmi = deg;
    for i in 1..=deg {
        psav[i] *= coef;
        dmi -= 1;
        coef *= dmi as f64;
    }
    // OCCT BuildCache non-rational branch (L1581-1589): scale by
    // LocalValue accumulated as LocalValue *= SpanDomain/ii (= 1/ii),
    // giving CachePoles[ii] = psav[ii-1] / (ii-1)!.
    let mut lv = 1.0f64;
    let mut out = Vec::with_capacity(deg + 1);
    for ii in 0..=deg {
        out.push(psav[ii] * lv);
        lv *= 1.0 / (ii as f64 + 1.0);
    }
    out
}

/// 2D variant of bezier_poles_to_coeffs (BSplCLib_BzSyntaxes.cxx L75-85).
fn bezier_poles_to_coeffs_2d(poles: &[DVec2]) -> Vec<DVec2> {
    let p3: Vec<DVec3> = poles.iter().map(|p| DVec3::new(p.x, p.y, 0.0)).collect();
    let c3 = bezier_poles_to_coeffs(&p3);
    c3.iter().map(|c| DVec2::new(c.x, c.y)).collect()
}

/// OCCT BSplCLib::IncreaseDegree for a Bezier pole sequence: raise the
/// degree of a Bezier curve by inserting the pole combination formula.
pub fn bezier_increase_degree(poles: &[DVec3], new_deg: usize) -> Vec<DVec3> {
    let mut p = poles.to_vec();
    while p.len() - 1 < new_deg {
        let deg = p.len() - 1;
        let mut np = Vec::with_capacity(p.len() + 1);
        np.push(p[0]);
        for i in 1..p.len() {
            let t = i as f64 / (deg + 1) as f64;
            np.push(p[i - 1] * t + p[i] * (1.0 - t));
        }
        np.push(p[p.len() - 1]);
        p = np;
    }
    p
}

/// OCCT BSplCLib::IncreaseDegree for 2D Bezier pole sequences.
pub fn bezier_increase_degree_2d(poles: &[DVec2], new_deg: usize) -> Vec<DVec2> {
    let mut p = poles.to_vec();
    while p.len() - 1 < new_deg {
        let deg = p.len() - 1;
        let mut np = Vec::with_capacity(p.len() + 1);
        np.push(p[0]);
        for i in 1..p.len() {
            let t = i as f64 / (deg + 1) as f64;
            np.push(p[i - 1] * t + p[i] * (1.0 - t));
        }
        np.push(p[p.len() - 1]);
        p = np;
    }
    p
}

/// OCCT GeomInt_TheMultiLineOfWLApprox + TheMultiLineToolOfWLApprox —
/// access to the WLine points of one approximation part (indices are
/// 1-based into the WLine point array, OCCT IntPatch_WLine 1-based).
pub struct WLineAccess<'a> {
    pub line: &'a IntPatchLine,
    pub indicemin: usize,
    pub indicemax: usize,
    pub nbp3d: usize,
    pub nbp2d: usize,
    pub approx_u1v1: bool,
    pub approx_u2v2: bool,
    pub p2d_on_first: bool,
    pub xo: f64,
    pub yo: f64,
    pub zo: f64,
    pub u1o: f64,
    pub v1o: f64,
    pub u2o: f64,
    pub v2o: f64,
    // The two surfaces of the face pair — used by the SvSurfaces tangency
    // (OCCT GeomInt_WLApprox with the quadric implicit surface).
    pub s1: &'a Surface3,
    pub s2: &'a Surface3,
    // The UV domains of the two surfaces (the PSurf domain used by
    // FillInitialVectorOfSolution, ApproxInt_ImpPrmSvSurfaces.gxx L865-1039).
    pub uv1: [f64; 4],
    pub uv2: [f64; 4],
}

impl<'a> WLineAccess<'a> {
    pub fn first_point(&self) -> usize {
        self.indicemin
    }
    pub fn last_point(&self) -> usize {
        self.indicemax
    }
    pub fn nb_p3d(&self) -> usize {
        self.nbp3d
    }
    pub fn nb_p2d(&self) -> usize {
        self.nbp2d
    }
    /// OCCT ApproxInt_MultiLine::WhatStatus (ApproxInt_MultiLine.gxx
    /// L154-160): a non-NULL SvSurfaces pointer means extra points can be
    /// added.  ApproxInt_Approx::buildCurve always passes the SvSurfaces
    /// pointer (ApproxInt_Approx.gxx L638-654), so the status is PointsAdded.
    pub fn what_status(&self) -> ApproxStatus {
        ApproxStatus::PointsAdded
    }
    fn wp(&self, index: usize) -> &WLinePnt {
        &self.line.wline_pnts[index - 1]
    }
    /// OCCT MultiLine::Value(Index, TabPnt): the 3D point (translated).
    pub fn value_p3d(&self, index: usize) -> DVec3 {
        let p = self.wp(index);
        DVec3::new(p.p3d.x + self.xo, p.p3d.y + self.yo, p.p3d.z + self.zo)
    }
    /// OCCT MultiLine::Value(Index, TabPnt2d): the 2D points (translated).
    pub fn value_p2d(&self, index: usize) -> Vec<DVec2> {
        let p = self.wp(index);
        let mut out = Vec::with_capacity(self.nbp2d);
        if self.nbp2d == 1 {
            if self.p2d_on_first {
                out.push(DVec2::new(p.u1 + self.u1o, p.v1 + self.v1o));
            } else {
                out.push(DVec2::new(p.u2 + self.u2o, p.v2 + self.v2o));
            }
        } else {
            out.push(DVec2::new(p.u1 + self.u1o, p.v1 + self.v1o));
            if self.nbp2d >= 2 {
                out.push(DVec2::new(p.u2 + self.u2o, p.v2 + self.v2o));
            }
        }
        out
    }
    /// OCCT GeomInt_TheMultiLineOfWLApprox::Tangency (ApproxInt_MultiLine.gxx
    /// L224-298) with the SvSurfaces.  For the quadric case
    /// (ApproxInt_ImpPrmSvSurfaces::Compute, gxx L437-765): the implicit
    /// quadric is the first quadric surface (Perform L246-306: typeS1 quadric
    /// -> Quad = S1, PSurf = S2), the intersection tangent is
    /// Tg = N_imp x N_prm (L693) normalized (L706), and the 2D tangents come
    /// from NonSingularProcessing (L287-321).  The singular cases return None
    /// so the constraint degrades to PassPoint (rcad LeastSquare::affect).
    pub fn tangency(&self, index: usize) -> Option<(Vec<DVec3>, Vec<DVec2>)> {
        let p = self.wp(index);
        let is_s1_quad = matches!(
            self.s1,
            Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Cone(_)
        );
        let is_s2_quad = matches!(
            self.s2,
            Surface3::Plane(_) | Surface3::Cylinder(_) | Surface3::Sphere(_) | Surface3::Cone(_)
        );
        if !is_s1_quad && !is_s2_quad {
            // No quadric — the parametric-parametric SvSurfaces path is not
            // ported (the tangency is unavailable).
            return None;
        }
        let (imp_surf, imp_u, imp_v) = if is_s1_quad {
            (self.s1, p.u1, p.v1)
        } else {
            (self.s2, p.u2, p.v2)
        };
        let (prm_surf, prm_u, prm_v) = if is_s1_quad {
            (self.s2, p.u2, p.v2)
        } else {
            (self.s1, p.u1, p.v1)
        };
        // OCCT ApproxInt_ImpPrmSvSurfaces::FillInitialVectorOfSolution
        // (L865-1039): the PSurf parameters must lie inside the PSurf domain
        // (with the 1e-10 slack); an out-of-domain non-periodic axis rejects the
        // point so the tangency constraint degrades to PassPoint.  The PSurf is
        // the parametric surface (the second one when the quadric is first).
        let prm_uv = if is_s1_quad { self.uv2 } else { self.uv1 };
        let (prm_binfu, prm_bsupu, prm_binfv, prm_bsupv) =
            (prm_uv[0], prm_uv[1], prm_uv[2], prm_uv[3]);
        let prm_u_per = if prm_surf.is_u_periodic() {
            Some(2.0 * std::f64::consts::PI)
        } else {
            None
        };
        let prm_v_per = if prm_surf.is_v_periodic() {
            Some(2.0 * std::f64::consts::PI)
        } else {
            None
        };
        // OCCT L886-937: u out of [binfu-1e-10, bsupu+1e-10] wraps only when
        // periodic, otherwise the initial solution is rejected.
        let mut trans_u = 0.0;
        if prm_u < prm_binfu - 0.0000000001 {
            match prm_u_per {
                Some(d) => {
                    while prm_u + trans_u < prm_binfu {
                        trans_u += d;
                    }
                }
                None => return None,
            }
        } else if prm_u > prm_bsupu + 0.0000000001 {
            match prm_u_per {
                Some(d) => {
                    while prm_u + trans_u > prm_bsupu {
                        trans_u -= d;
                    }
                }
                None => return None,
            }
        }
        let mut trans_v = 0.0;
        if prm_v < prm_binfv - 0.0000000001 {
            match prm_v_per {
                Some(d) => {
                    while prm_v + trans_v < prm_binfv {
                        trans_v += d;
                    }
                }
                None => return None,
            }
        } else if prm_v > prm_bsupv + 0.0000000001 {
            match prm_v_per {
                Some(d) => {
                    while prm_v + trans_v > prm_bsupv {
                        trans_v -= d;
                    }
                }
                None => return None,
            }
        }
        // OCCT L938-939/L1015-1016: X = params + Translation.
        let prm_ux = prm_u + trans_u;
        let prm_vx = prm_v + trans_v;
        // The implicit quadric's normal at the point (aQSurf.Normale(MyPnt),
        // ImpPrmSvSurfaces.gxx L625) and the parametric surface's normal
        // (ThePSurfaceTool::D1 L591).
        let n_imp = imp_surf.normal_at(imp_u, imp_v).normalize_or_zero();
        let (_, du_prm, dv_prm) = prm_surf.derivatives(prm_ux, prm_vx);
        let n_prm = du_prm.cross(dv_prm).normalize_or_zero();
        if n_imp.length_squared() < 1.0e-14 || n_prm.length_squared() < 1.0e-14 {
            return None;
        }
        // Tg = N_imp x N_prm (L693), normalized (L706).
        let tg = n_imp.cross(n_prm).normalize_or_zero();
        if tg.length_squared() < 1.0e-14 {
            return None;
        }
        // The 2D tangents on both surfaces (NonSingularProcessing L287-321;
        // the singular case degrades to PassPoint).
        let (_, du_imp, dv_imp) = imp_surf.derivatives(imp_u, imp_v);
        let normal_imp = du_imp.cross(dv_imp);
        if normal_imp.length_squared() < 1.0e-14 {
            return None;
        }
        let normal_prm = du_prm.cross(dv_prm);
        if normal_prm.length_squared() < 1.0e-14 {
            return None;
        }
        let (ts_imp, ts_prm) = (
            nonsingular_tangent_2d(du_imp, dv_imp, normal_imp, tg),
            nonsingular_tangent_2d(du_prm, dv_prm, normal_prm, tg),
        );
        let (ts1, ts2) = if is_s1_quad { (ts_imp, ts_prm) } else { (ts_prm, ts_imp) };
        if std::env::var("RCAD_TG_DEBUG").is_ok() {
            eprintln!("[TG] idx={} p=({:.9},{:.9},{:.9}) tg=({:.9},{:.9},{:.9}) ts1=({:.6},{:.6}) ts2=({:.6},{:.6}) s1={:?} s2={:?} u1={:.6} v1={:.6} u2={:.6} v2={:.6} is_s1_quad={}", index, p.p3d.x, p.p3d.y, p.p3d.z, tg.x, tg.y, tg.z, ts1.x, ts1.y, ts2.x, ts2.y, surface_kind(self.s1), surface_kind(self.s2), p.u1, p.v1, p.u2, p.v2, is_s1_quad);
        }
        Some((vec![tg], vec![ts1, ts2]))
    }
    /// OCCT LineTool::MakeMLBetween / MakeMLOneMorePoint — no SvSurfaces ->
    /// no point insertion possible.
    pub fn make_ml_one_more_point(
        &self,
        _low: usize,
        _high: usize,
        _indbad: usize,
    ) -> Option<WLineAccess<'a>> {
        None
    }
    pub fn make_ml_between(&self, _low: usize, _high: usize, _n: usize) -> Option<WLineAccess<'a>> {
        None
    }
}

/// OCCT NonSingularProcessing (ApproxInt_ImpPrmSvSurfaces.gxx L287-321): the
/// 2D tangent (theTg2D) on a surface with the derivative basis (DU, DV) such
/// that Tg3D = DU*Tg2D.X() + DV*Tg2D.Y() holds.  aNormal = DU x DV must be
/// non-zero (checked by the caller).
fn nonsingular_tangent_2d(du: DVec3, dv: DVec3, normal: DVec3, tg3d: DVec3) -> DVec2 {
    // If T = A*U + B*V then
    //   A x T = (A x B)*V
    //   B x T = (B x A)*U
    let tgu = tg3d.cross(du);
    let tgv = tg3d.cross(dv);
    let sq_magn = normal.length_squared();
    let delta_u = tgv.length_squared() / sq_magn;
    let delta_v = tgu.length_squared() / sq_magn;
    DVec2::new(
        delta_u.sqrt().copysign(tgv.dot(normal)),
        -delta_v.sqrt().copysign(tgu.dot(normal)),
    )
}

/// Debug helper: a short tag for the surface kind.
fn surface_kind(s: &Surface3) -> &'static str {
    match s {
        Surface3::Plane(_) => "Plane",
        Surface3::Cylinder(_) => "Cylinder",
        Surface3::Sphere(_) => "Sphere",
        Surface3::Cone(_) => "Cone",
        Surface3::Torus(_) => "Torus",
        _ => "Other",
    }
}

/// OCCT Approx_ComputeLine::CheckMultiCurve (Approx_ComputeLine.gxx
/// L134-428): reject a fit whose poles make a loop that the data does not
/// justify; on rejection `the_indbad` locates the longest segment.
fn check_multicurve(
    the_multi_curve: &MultiCurve,
    ml: &WLineAccess,
    the_indfirst: usize,
    the_indlast: usize,
    the_indbad: &mut usize,
) -> bool {
    let nbp3d = ml.nb_p3d();
    let nbp2d = ml.nb_p2d();
    let coeff = 4.0;
    if nbp3d > 1 {
        return true;
    }
    let min_scal_prod = -0.9;
    let sq_tol3d = 1.0e-14; // Precision::SquareConfusion
    *the_indbad = 0;
    let mut indbads = [0usize; 4];
    let nb_cur = the_multi_curve.nb_curves();
    let mut loop_found = false;
    if the_multi_curve.dimension(1) == 3 {
        let mut a_poles = Vec::new();
        the_multi_curve.curve(1, &mut a_poles);
        let mut first_vec = DVec3::ZERO;
        let mut indp = 2usize;
        while indp <= a_poles.len() {
            first_vec = a_poles[indp - 1] - a_poles[0];
            indp += 1;
            let a_length = first_vec.length();
            if a_length > 1.0e-12 {
                first_vec /= a_length;
                break;
            }
        }
        let mut mid_pnt = a_poles[indp - 1];
        while indp <= a_poles.len() {
            let mut second_vec = a_poles[indp - 1] - mid_pnt;
            let a_length = second_vec.length();
            if a_length <= 1.0e-12 {
                indp += 1;
                continue;
            }
            second_vec /= a_length;
            let scal_prod = first_vec.dot(second_vec);
            if scal_prod < min_scal_prod {
                loop_found = true;
                break;
            }
            first_vec = second_vec;
            mid_pnt = a_poles[indp - 1];
            indp += 1;
        }
        if loop_found {
            for first_ind in the_indfirst..=the_indlast - 2 {
                let first_pnt = ml.value_p3d(first_ind);
                let mut real_loop = true;
                for k in first_ind + 1..the_indlast {
                    let pnt1 = ml.value_p3d(k);
                    let pnt2 = ml.value_p3d(k + 1);
                    if first_pnt.distance_squared(pnt1) <= sq_tol3d
                        || first_pnt.distance_squared(pnt2) <= sq_tol3d
                    {
                        loop_found = false;
                        real_loop = false;
                        break;
                    }
                    let vec1 = (pnt1 - first_pnt).normalize_or_zero();
                    let vec2 = (pnt2 - first_pnt).normalize_or_zero();
                    if vec1.dot(vec2) < min_scal_prod {
                        loop_found = false;
                        real_loop = false;
                        break;
                    }
                }
                if !real_loop {
                    break;
                }
            }
        }
        if loop_found {
            let mut max_sq_dist = 0.0f64;
            let mut min_sq_dist = f64::INFINITY;
            for k in the_indfirst + 1..=the_indlast {
                let prev_pnt = ml.value_p3d(k - 1);
                let cur_pnt = ml.value_p3d(k);
                let a_sq_dist = prev_pnt.distance_squared(cur_pnt);
                if a_sq_dist > max_sq_dist {
                    max_sq_dist = a_sq_dist;
                    indbads[1] = k;
                }
                if a_sq_dist > 1.0e-12 && a_sq_dist < min_sq_dist {
                    min_sq_dist = a_sq_dist;
                }
            }
            let relation = max_sq_dist / min_sq_dist;
            if relation < coeff {
                loop_found = false;
            } else {
                for indcur in 2..=nb_cur {
                    max_sq_dist = 0.0;
                    for k in the_indfirst + 1..=the_indlast {
                        let prev_pnt = ml.value_p2d(k - 1);
                        let cur_pnt = ml.value_p2d(k);
                        let a_sq_dist = (prev_pnt[indcur - 1] - cur_pnt[indcur - 1]).length_squared();
                        if a_sq_dist > max_sq_dist {
                            max_sq_dist = a_sq_dist;
                            indbads[indcur] = k;
                        }
                    }
                }
            }
        }
    } else {
        // 2d case.
        let mut a_poles2d = Vec::new();
        the_multi_curve.curve2d(1, &mut a_poles2d);
        let a_sq_norm_toler = 2.220446049250313e-16 * 2.220446049250313e-16;
        let mut first_vec = a_poles2d[1] - a_poles2d[0];
        let mut a_vec_sq_norm = first_vec.length_squared();
        if a_vec_sq_norm < a_sq_norm_toler {
            *the_indbad = the_indfirst + 1;
            return false;
        }
        first_vec /= a_sq_norm_toler.sqrt();
        let mut mid_pnt = a_poles2d[1];
        for k in 2..a_poles2d.len() {
            let second_vec = a_poles2d[k] - mid_pnt;
            a_vec_sq_norm = second_vec.length_squared();
            if a_vec_sq_norm < a_sq_norm_toler {
                *the_indbad = the_indfirst + k - 1 + 1;
                return false;
            }
            let second_vec = second_vec / a_vec_sq_norm.sqrt();
            let scal_prod = first_vec.dot(second_vec);
            if scal_prod < min_scal_prod {
                loop_found = true;
                break;
            }
            first_vec = second_vec;
            mid_pnt = a_poles2d[k];
        }
        if loop_found {
            for first_ind in the_indfirst..=the_indlast - 2 {
                let first_pnt = ml.value_p2d(first_ind)[0];
                let mut real_loop = true;
                for k in first_ind + 1..the_indlast {
                    let pnt1 = ml.value_p2d(k)[0];
                    let pnt2 = ml.value_p2d(k + 1)[0];
                    if first_pnt.distance_squared(pnt1) <= sq_tol3d
                        || first_pnt.distance_squared(pnt2) <= sq_tol3d
                    {
                        loop_found = false;
                        real_loop = false;
                        break;
                    }
                    let vec1 = (pnt1 - first_pnt).normalize_or_zero();
                    let vec2 = (pnt2 - first_pnt).normalize_or_zero();
                    if vec1.dot(vec2) < min_scal_prod {
                        loop_found = false;
                        real_loop = false;
                        break;
                    }
                }
                if !real_loop {
                    break;
                }
            }
        }
        if loop_found {
            for indcur in 1..=nb_cur {
                let mut max_sq_dist = 0.0f64;
                let mut min_sq_dist = f64::INFINITY;
                for k in the_indfirst + 1..=the_indlast {
                    let prev_pnt = ml.value_p2d(k - 1);
                    let cur_pnt = ml.value_p2d(k);
                    let a_sq_dist = (prev_pnt[indcur - 1] - cur_pnt[indcur - 1]).length_squared();
                    if a_sq_dist > max_sq_dist {
                        max_sq_dist = a_sq_dist;
                        indbads[indcur] = k;
                    }
                    if a_sq_dist > 1.0e-12 && a_sq_dist < min_sq_dist {
                        min_sq_dist = a_sq_dist;
                    }
                }
                let relation = max_sq_dist / min_sq_dist;
                if relation < coeff {
                    loop_found = false;
                }
            }
        }
    }
    for i in 1..=3 {
        if indbads[i] != 0 {
            *the_indbad = indbads[i];
            break;
        }
    }
    if !loop_found {
        *the_indbad = 0;
    }
    !loop_found
}

// ============================================================================
// Approx_MCurvesToBSpCurve (Approx_MCurvesToBSpCurve.cxx L58-292) +
// Convert_CompBezierCurvesToBSplineCurveBase (Convert_...Base.hxx L50-173) —
// concatenate the per-interval Bezier MultiCurves into one BSpline.
// ============================================================================

/// OCCT Convert_CompBezierCurvesToBSplineCurveBase::Perform — merge a
/// sequence of adjacent Bezier pole tables into one BSpline (degree, knots,
/// multiplicities, poles).  Returns (degree, knots, mults, poles).
pub fn convert_comp_bezier(poles_seq: &[Vec<DVec3>]) -> (usize, Vec<f64>, Vec<usize>, Vec<DVec3>) {
    let mut curve_poles: Vec<DVec3> = Vec::new();
    let mut curve_knots: Vec<f64> = Vec::new();
    let mut knots_mults: Vec<usize> = Vec::new();
    let nb_curv = poles_seq.len();
    let mut degree = 0usize;
    for s in poles_seq {
        degree = degree.max(s.len() - 1);
    }
    let mut kn_vals = vec![0.0f64; nb_curv];
    let mut a_p1 = DVec3::ZERO;
    let mut a_det = 0.0;
    for (i, seg) in poles_seq.iter().enumerate() {
        // 1- raise the Bezier curve to the maximum degree.
        let points = bezier_increase_degree(seg, degree);
        // 2- process the junction node.
        if i == 0 {
            for j in 0..degree {
                curve_poles.push(points[j]);
            }
            kn_vals[0] = 1.0;
            knots_mults.push(degree + 1);
            a_det = 1.0;
        }
        if i != 0 {
            let a_p2 = points[0];
            let a_p3 = points[1];
            let a_v1 = a_p2 - a_p1;
            let a_v2 = a_p3 - a_p2;
            let a_d1 = a_v1.length_squared();
            let a_d2 = a_v2.length_squared();
            let angular = 1.0e-4;
            if degree > 1 && a_d1 > 1.0e-12 && a_d2 > 1.0e-12 && v_parallel(&a_v1, &a_v2, angular) {
                let a_lambda = (a_d2 / a_d1).sqrt();
                if kn_vals[i - 1] * a_lambda > 10.0 * eps_val(a_det) {
                    knots_mults.push(degree - 1);
                    kn_vals[i] = kn_vals[i - 1] * a_lambda;
                } else {
                    curve_poles.push(points[0]);
                    knots_mults.push(degree);
                    kn_vals[i] = 1.0;
                }
            } else {
                curve_poles.push(points[0]);
                knots_mults.push(degree);
                kn_vals[i] = 1.0;
            }
            a_det += kn_vals[i];
            for j in 1..degree {
                curve_poles.push(points[j]);
            }
        }
        if i == nb_curv - 1 {
            curve_poles.push(points[degree]);
            knots_mults.push(degree + 1);
        }
        a_p1 = points[degree];
    }
    // Correct the nodal values to [0, 1].
    curve_knots.push(0.0);
    for i in 1..nb_curv {
        curve_knots.push(curve_knots[i - 1] + kn_vals[i - 1] / a_det);
    }
    curve_knots.push(1.0);
    (degree, curve_knots, knots_mults, curve_poles)
}

fn v_parallel(v1: &DVec3, v2: &DVec3, angular: f64) -> bool {
    let c = v1.cross(*v2);
    c.length() <= angular * v1.length() * v2.length()
}

fn eps_val(x: f64) -> f64 {
    x.abs() * 2.220446049250313e-16
}

/// OCCT Approx_MCurvesToBSpCurve.
pub struct MCurvesToBSpCurve {
    pub my_curves: Vec<MultiCurve>,
    pub my_spline: MultiBSpCurve,
    pub my_done: bool,
}

impl MCurvesToBSpCurve {
    pub fn new() -> Self {
        MCurvesToBSpCurve {
            my_curves: Vec::new(),
            my_spline: MultiBSpCurve {
                degree: 0,
                knots: Vec::new(),
                mults: Vec::new(),
                poles: Vec::new(),
            },
            my_done: false,
        }
    }
    pub fn reset(&mut self) {
        self.my_done = false;
        self.my_curves.clear();
    }
    pub fn append(&mut self, mc: MultiCurve) {
        self.my_curves.push(mc);
    }
    pub fn perform(&mut self) {
        let the_seq = &self.my_curves;
        let nbcu = the_seq.len();
        if nbcu == 1 {
            let cu = &the_seq[0];
            let deg = cu.degree();
            let knots = vec![0.0, 1.0];
            let mults = vec![deg + 1, deg + 1];
            self.my_spline = MultiBSpCurve::from_bezier(cu, knots, mults);
        } else {
            let p = the_seq[nbcu - 1].value(1);
            let nb3d = p.nb_points();
            let nb2d = p.nb_points2d();
            let mut the_poles_spl: Vec<DVec3> = Vec::new();
            let mut the_poles_spl2d: Vec<DVec2> = Vec::new();
            let mut the_knots: Vec<f64> = Vec::new();
            let mut the_mults: Vec<usize> = Vec::new();
            let mut deg = 0usize;
            if nb3d != 0 {
                let mut seq: Vec<Vec<DVec3>> = Vec::new();
                for cu in the_seq {
                    let mut p3 = Vec::new();
                    cu.curve(1, &mut p3);
                    seq.push(p3);
                }
                let (d, k, m, pl) = convert_comp_bezier(&seq);
                the_knots = k;
                the_mults = m;
                the_poles_spl = pl;
                deg = d;
            } else if nb2d != 0 {
                // 2D-only: convert the first 2D curve.
                let mut seq: Vec<Vec<DVec3>> = Vec::new();
                for cu in the_seq {
                    let mut p3 = Vec::new();
                    cu.curve2d(1, &mut p3);
                    seq.push(p3.iter().map(|q| DVec3::new(q.x, q.y, 0.0)).collect());
                }
                let (d, k, m, pl) = convert_comp_bezier(&seq);
                the_knots = k;
                the_mults = m;
                the_poles_spl = pl;
                deg = d;
                the_poles_spl2d = the_poles_spl.iter().map(|q| DVec2::new(q.x, q.y)).collect();
            }
            let nb_poles_spl = if nb3d != 0 { the_poles_spl.len() } else { the_poles_spl2d.len() };
            let mut tab_mu: Vec<MultiPoint> = (0..nb_poles_spl)
                .map(|_| MultiPoint::new(nb3d, nb2d))
                .collect();
            if nb3d != 0 {
                for (j, mp) in tab_mu.iter_mut().enumerate() {
                    mp.set_point(1, the_poles_spl[j]);
                }
            } else if nb2d != 0 {
                for (j, mp) in tab_mu.iter_mut().enumerate() {
                    mp.set_point2d(1 + nb3d, the_poles_spl2d[j]);
                }
            }
            let thefirst = if nb3d != 0 { 1 } else { 2 };
            let mut kpoles3d = 0usize;
            let mut kpoles2d = 0usize;
            let mut kpol = 0usize;
            for (i, cu) in the_seq.iter().enumerate() {
                let mydegre = cu.degree();
                let last = if the_mults[i + 1] == deg {
                    deg + 1 // C0
                } else {
                    deg // C1
                };
                let last = if i == nbcu - 1 { deg + 1 } else { last };
                let first = if i == 0 {
                    1
                } else if the_mults[i] == deg - 1 || the_mults[i] == deg {
                    2
                } else {
                    1
                };
                // 3D curves 2..nb3d (unused when nb3d == 1).
                for j in 2..=nb3d {
                    kpol = kpoles3d;
                    let mut the_poles = Vec::new();
                    cu.curve(j, &mut the_poles);
                    let inc = deg - mydegre;
                    let points = if inc > 0 {
                        bezier_increase_degree(&the_poles, deg)
                    } else {
                        the_poles
                    };
                    for k in first..=last {
                        tab_mu[kpol].set_point(j, points[k - 1]);
                        kpol += 1;
                    }
                }
                kpoles3d = kpol;
                // 2D curves thefirst..nb2d.
                for j in thefirst..=nb2d {
                    kpol = kpoles2d;
                    let mut the_poles2d = Vec::new();
                    cu.curve2d(j + nb3d, &mut the_poles2d);
                    let inc = deg - mydegre;
                    let points2d = if inc > 0 {
                        bezier_increase_degree_2d(&the_poles2d, deg)
                    } else {
                        the_poles2d
                    };
                    for k in first..=last {
                        tab_mu[kpol].set_point2d(j + nb3d, points2d[k - 1]);
                        kpol += 1;
                    }
                }
                kpoles2d = kpol;
            }
            self.my_spline = MultiBSpCurve::from_multipoints(tab_mu, the_knots, the_mults, deg);
        }
        self.my_done = true;
    }
    pub fn value(&self) -> MultiBSpCurve {
        self.my_spline.clone()
    }
}

// ============================================================================
// GeomInt_WLApprox / ApproxInt_Approx (GeomInt_WLApprox.hxx + _0.cxx +
// ApproxInt_Approx.gxx L166-751) — the top of the WLine approximation.
// ============================================================================

/// OCCT GeomInt_WLApprox (the ApproxInt_Approx template instantiated with
/// the WLine multi-line).
pub struct WLineApprox {
    // Approx_Data.
    pub my_bezier_approx: bool,
    pub xo: f64,
    pub yo: f64,
    pub zo: f64,
    pub u1o: f64,
    pub v1o: f64,
    pub u2o: f64,
    pub v2o: f64,
    pub approx_xyz: bool,
    pub approx_u1v1: bool,
    pub approx_u2v2: bool,
    pub indicemin: usize,
    pub indicemax: usize,
    pub my_nb_pnt_max: usize,
    pub parametrization: ApproxParamType,
    // ApproxInt_Approx members.
    pub my_compute_line: ComputeLine,
    pub my_compute_line_bezier: ComputeLine,
    pub my_bez_to_bspl: MCurvesToBSpCurve,
    pub my_with_tangency: bool,
    pub my_tol3d: f64,
    pub my_tol2d: f64,
    pub my_deg_min: usize,
    pub my_deg_max: usize,
    pub my_nb_iter_max: i32,
    pub my_tol_reached3d: f64,
    pub my_tol_reached2d: f64,
    pub my_knots: Vec<usize>,
    pub my_done: bool,
    pub my_value: MultiBSpCurve,
}

impl WLineApprox {
    pub fn new() -> Self {
        WLineApprox {
            my_bezier_approx: true,
            xo: 0.0,
            yo: 0.0,
            zo: 0.0,
            u1o: 0.0,
            v1o: 0.0,
            u2o: 0.0,
            v2o: 0.0,
            approx_xyz: true,
            approx_u1v1: true,
            approx_u2v2: true,
            indicemin: 0,
            indicemax: 0,
            my_nb_pnt_max: 30,
            parametrization: ApproxParamType::ChordLength,
            my_compute_line: ComputeLine::new(4, 8, 0.001, 0.001, 5, true, ApproxParamType::ChordLength, false),
            my_compute_line_bezier: ComputeLine::new(4, 8, 0.001, 0.001, 5, true, ApproxParamType::ChordLength, false),
            my_bez_to_bspl: MCurvesToBSpCurve::new(),
            my_with_tangency: true,
            my_tol3d: 0.001,
            my_tol2d: 0.001,
            my_deg_min: 4,
            my_deg_max: 8,
            my_nb_iter_max: 5,
            my_tol_reached3d: 0.0,
            my_tol_reached2d: 0.0,
            my_knots: Vec::new(),
            my_done: false,
            my_value: MultiBSpCurve {
                degree: 0,
                knots: Vec::new(),
                mults: Vec::new(),
                poles: Vec::new(),
            },
        }
    }

    /// OCCT ApproxInt_Approx::SetParameters (L399-427).
    pub fn set_parameters(
        &mut self,
        tol3d: f64,
        tol2d: f64,
        deg_min: usize,
        deg_max: usize,
        nb_iter_max: i32,
        nb_pnt_max: usize,
        approx_with_tangency: bool,
        parametrization: ApproxParamType,
    ) {
        let ratio_tol = 1.5;
        self.my_nb_pnt_max = nb_pnt_max;
        self.my_with_tangency = approx_with_tangency;
        self.my_tol3d = tol3d / ratio_tol;
        self.my_tol2d = tol2d / ratio_tol;
        self.my_deg_min = deg_min;
        self.my_deg_max = deg_max;
        self.my_nb_iter_max = nb_iter_max;
        self.my_compute_line.init(deg_min, deg_max, self.my_tol3d, self.my_tol2d, nb_iter_max, true, parametrization);
        self.my_compute_line_bezier.init(deg_min, deg_max, self.my_tol3d, self.my_tol2d, nb_iter_max, true, parametrization);
        if !approx_with_tangency {
            self.my_compute_line.set_constraints(AppParConstraint::PassPoint, AppParConstraint::PassPoint);
            self.my_compute_line_bezier.set_constraints(AppParConstraint::PassPoint, AppParConstraint::PassPoint);
        }
        self.my_bezier_approx = true;
    }

    /// OCCT ApproxInt_Approx::Perform (L184-218) — the shared core.
    pub fn perform(
        &mut self,
        ml: &WLineAccess,
        approx_xyz: bool,
        approx_u1v1: bool,
        approx_u2v2: bool,
        indicemin: usize,
        indicemax: usize,
    ) {
        // prepareDS (L521-534).
        self.my_tol_reached3d = 0.0;
        self.my_tol_reached2d = 0.0;
        self.approx_u1v1 = approx_u1v1;
        self.approx_u2v2 = approx_u2v2;
        self.approx_xyz = approx_xyz;
        self.indicemin = indicemin;
        self.indicemax = indicemax;
        self.parametrization = self.my_compute_line_bezier.par;
        let nbpntbez = self.indicemax - self.indicemin;
        self.my_bezier_approx = nbpntbez >= 5; // aMinNbPointsForApprox
        // fillData (L501-517).
        self.fill_data(ml);
        // buildKnots (L538-619).
        self.build_knots(ml);
        if self.my_knots.len() == 2 && self.indicemax - self.indicemin > 2 * self.my_nb_pnt_max {
            self.my_knots[1] = (self.indicemax - self.indicemin) / 2;
            self.my_knots.push(self.indicemax);
        }
        self.my_compute_line_bezier.init(
            self.my_deg_min,
            self.my_deg_max,
            self.my_tol3d,
            self.my_tol2d,
            self.my_nb_iter_max,
            true,
            self.parametrization,
        );
        self.build_curve(ml);
        self.my_done = self.my_compute_line_bezier.nb_multi_curves() > 0;
    }

    /// OCCT ApproxInt_Approx::fillData (L501-517) — ComputeTrsf3d/2d.
    /// The translation minima are taken over the WHOLE WLine (ComputeTrsf3d
    /// L35-53 and ComputeTrsf2d L57-85 iterate theLine->NbPnts()).
    fn fill_data(&mut self, ml: &WLineAccess) {
        let all_pnts: Vec<WLinePnt> = ml.line.wline_pnts.clone();
        // ComputeTrsf3d.
        if self.approx_xyz {
            let mut xmin = f64::INFINITY;
            let mut ymin = f64::INFINITY;
            let mut zmin = f64::INFINITY;
            for p in &all_pnts {
                xmin = xmin.min(p.p3d.x);
                ymin = ymin.min(p.p3d.y);
                zmin = zmin.min(p.p3d.z);
            }
            self.xo = -xmin;
            self.yo = -ymin;
            self.zo = -zmin;
        } else {
            self.xo = 0.0;
            self.yo = 0.0;
            self.zo = 0.0;
        }
        // ComputeTrsf2d on surface 1.
        if self.approx_u1v1 {
            let mut umin = f64::INFINITY;
            let mut vmin = f64::INFINITY;
            for p in &all_pnts {
                umin = umin.min(p.u1);
                vmin = vmin.min(p.v1);
            }
            self.u1o = -umin;
            self.v1o = -vmin;
        } else {
            self.u1o = 0.0;
            self.v1o = 0.0;
        }
        // ComputeTrsf2d on surface 2.
        if self.approx_u2v2 {
            let mut umin = f64::INFINITY;
            let mut vmin = f64::INFINITY;
            for p in &all_pnts {
                umin = umin.min(p.u2);
                vmin = vmin.min(p.v2);
            }
            self.u2o = -umin;
            self.v2o = -vmin;
        } else {
            self.u2o = 0.0;
            self.v2o = 0.0;
        }
    }

    /// OCCT ApproxInt_Approx::buildKnots (L538-619).
    fn build_knots(&mut self, ml: &WLineAccess) {
        self.my_knots.clear();
        if !self.my_bezier_approx {
            self.my_knots.push(self.indicemin);
            self.my_knots.push(self.indicemax);
            return;
        }
        let nbp3d = if self.approx_xyz { 1 } else { 0 };
        let nbp2d = (if self.approx_u1v1 { 1 } else { 0 }) + (if self.approx_u2v2 { 1 } else { 0 });
        let a_dim = nbp3d * 3 + nbp2d * 2;
        // Collect the part points into the concatenated coordinate array.
        let n = self.indicemax - self.indicemin + 1;
        let mut a_coords = vec![0.0f64; n * a_dim];
        for (i, idx) in (self.indicemin..=self.indicemax).enumerate() {
            let p3 = ml.value_p3d(idx);
            let p2 = ml.value_p2d(idx);
            let mut j = i * a_dim;
            if nbp3d > 0 {
                a_coords[j] = p3.x;
                a_coords[j + 1] = p3.y;
                a_coords[j + 2] = p3.z;
                j += 3;
            }
            if self.approx_u1v1 {
                a_coords[j] = p2[0].x;
                a_coords[j + 1] = p2[0].y;
                j += 2;
            }
            if self.approx_u2v2 {
                let k = if p2.len() >= 2 { 1 } else { 0 };
                a_coords[j] = p2[k].x;
                a_coords[j + 1] = p2[k].y;
                j += 2;
            }
        }
        let a_pars = approx_parameters(ml, self.indicemin, self.indicemax, self.parametrization);
        let mut knots = build_knots(&a_coords, a_dim, &a_pars.v, self.my_nb_pnt_max);
        // OCCT output knots are 1-based point indices; re-map to the actual
        // point indices (indicemin..).
        for k in knots.iter_mut() {
            *k += self.indicemin;
        }
        if knots.len() < 2 {
            knots = vec![self.indicemin, self.indicemax];
        }
        if std::env::var("RCAD_APX_DEBUG").is_ok() {
            eprintln!("  [KNOTS] indicemin={} indicemax={} knots={:?}", self.indicemin, self.indicemax, knots);
        }
        self.my_knots = knots;
    }

    /// OCCT ApproxInt_Approx::buildCurve (L623-751).
    fn build_curve(&mut self, ml: &WLineAccess) {
        self.my_bez_to_bspl.reset();
        let mut kind = 0usize;
        loop {
            let imin = self.my_knots[kind];
            let imax = self.my_knots[kind + 1];
            let sub = WLineAccess {
                line: ml.line,
                indicemin: imin,
                indicemax: imax,
                nbp3d: ml.nbp3d,
                nbp2d: ml.nbp2d,
                approx_u1v1: self.approx_u1v1,
                approx_u2v2: self.approx_u2v2,
                p2d_on_first: ml.p2d_on_first,
                xo: self.xo,
                yo: self.yo,
                zo: self.zo,
                u1o: self.u1o,
                v1o: self.v1o,
                u2o: self.u2o,
                v2o: self.v2o,
                s1: ml.s1,
                s2: ml.s2,
                uv1: ml.uv1,
                uv2: ml.uv2,
            };
            self.my_compute_line_bezier.perform(&sub);
            if self.my_compute_line_bezier.nb_multi_curves() == 0 {
                return;
            }
            self.update_tol_reached();
            // Transform the translated fits back (ApproxInt_Approx.gxx
            // L669-730: indice3d = 1, indice2d1 = 2, indice2d2 = 3, adjusted
            // by the approximated flags).
            let mut indice3d = 1usize;
            let mut indice2d1 = 2usize;
            let mut indice2d2 = 3usize;
            if !self.approx_xyz {
                indice2d1 -= 1;
                indice2d2 -= 1;
            }
            if !self.approx_u1v1 {
                indice2d2 -= 1;
            }
            let n = self.my_compute_line_bezier.nb_multi_curves();
            // OCCT ApproxInt_Approx.gxx L683-687/L699-706/L717-724: the
            // Translate/Transform pass runs in REVERSE order and modifies the
            // stored multi-curves IN PLACE (ChangeValue(...).Transform(...));
            // the Append pass (L735-738) then runs in FORWARD order on the
            // transformed curves.
            for idx in (0..n).rev() {
                let mc = self.my_compute_line_bezier.value_mut(idx + 1);
                if self.approx_xyz {
                    for mp in mc.poles.iter_mut() {
                        let p = &mut mp.p3d[indice3d - 1];
                        // OCCT L686: Transform(indice3d, -Xo, 1.0, -Yo, 1.0, -Zo, 1.0).
                        *p += DVec3::new(-self.xo, -self.yo, -self.zo);
                    }
                }
                if self.approx_u1v1 {
                    for mp in mc.poles.iter_mut() {
                        let loc = indice2d1 - 1 - self.approx_xyz as usize;
                        // OCCT L701: Transform2d(indice2d1, -U1o, 1.0, -V1o, 1.0).
                        mp.p2d[loc] += DVec2::new(-self.u1o, -self.v1o);
                    }
                }
                if self.approx_u2v2 {
                    for mp in mc.poles.iter_mut() {
                        let loc = indice2d2 - 1 - self.approx_xyz as usize;
                        if loc < mp.p2d.len() {
                            // OCCT L719: Transform2d(indice2d2, -U2o, 1.0, -V2o, 1.0).
                            mp.p2d[loc] += DVec2::new(-self.u2o, -self.v2o);
                        }
                    }
                }
            }
            // OCCT ApproxInt_Approx.gxx L735-738: the Append pass runs in
            // FORWARD order (1..NbMultiCurves).
            for idx in 1..=n {
                let mc = self.my_compute_line_bezier.value(idx).clone();
                self.my_bez_to_bspl.append(mc);
            }
            kind += 1;
            if kind < self.my_knots.len() - 1 {
                continue;
            }
            break;
        }
        self.my_bez_to_bspl.perform();
    }

    /// OCCT ApproxInt_Approx::UpdateTolReached (L438-455).
    fn update_tol_reached(&mut self) {
        let nb = self.my_compute_line_bezier.nb_multi_curves();
        for icur in 1..=nb {
            let (tol3d, tol2d) = self.my_compute_line_bezier.error(icur);
            self.my_tol_reached3d = self.my_tol_reached3d.max(tol3d);
            self.my_tol_reached2d = self.my_tol_reached2d.max(tol2d);
        }
    }

    pub fn tol_reached3d(&self) -> f64 {
        self.my_tol_reached3d * 1.5
    }
    pub fn tol_reached2d(&self) -> f64 {
        self.my_tol_reached2d * 1.5
    }
    pub fn is_done(&self) -> bool {
        self.my_done
    }
    pub fn nb_multi_curves(&self) -> usize {
        if self.my_bezier_approx {
            1
        } else {
            0
        }
    }
    pub fn value(&self) -> MultiBSpCurve {
        if self.my_bezier_approx {
            self.my_bez_to_bspl.value()
        } else {
            self.my_value.clone()
        }
    }
}

pub struct ComputeLine {
    pub mydegremin: usize,
    pub mydegremax: usize,
    pub mytol3d: f64,
    pub mytol2d: f64,
    pub par: ApproxParamType,
    pub mysquares: bool,
    pub mycut: bool,
    pub myitermax: i32,
    pub myfirstc: AppParConstraint,
    pub mylastc: AppParConstraint,
    pub myconstraints: Vec<ConstraintCouple>,
    pub myfirst_param: Option<VecD>,
    pub mymulti_curves: Vec<MultiCurve>,
    pub mypar: Vec<VecD>,
    pub tolers3d: Vec<f64>,
    pub tolers2d: Vec<f64>,
    pub the_multi_curve: MultiCurve,
    pub currenttol3d: f64,
    pub currenttol2d: f64,
    pub myparameters: VecD,
    pub alldone: bool,
    pub tolreached: bool,
    pub my_multi_line_nb: i32,
    pub my_is_clear: bool,
}

impl ComputeLine {
    /// OCCT Approx_ComputeLine constructor (L779-802) + Init (L1692-1709).
    pub fn new(
        degreemin: usize,
        degreemax: usize,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParamType,
        squares: bool,
    ) -> Self {
        ComputeLine {
            mydegremin: degreemin,
            mydegremax: degreemax,
            mytol3d: tolerance3d,
            mytol2d: tolerance2d,
            par: parametrization,
            mysquares: squares,
            mycut: cutting,
            myitermax: nb_iterations,
            myfirstc: AppParConstraint::TangencyPoint,
            mylastc: AppParConstraint::TangencyPoint,
            myconstraints: vec![
                ConstraintCouple { index: 0, constraint: AppParConstraint::TangencyPoint },
                ConstraintCouple { index: 0, constraint: AppParConstraint::TangencyPoint },
            ],
            myfirst_param: None,
            mymulti_curves: Vec::new(),
            mypar: Vec::new(),
            tolers3d: Vec::new(),
            tolers2d: Vec::new(),
            the_multi_curve: MultiCurve::new(2, 1, 2),
            currenttol3d: f64::INFINITY,
            currenttol2d: f64::INFINITY,
            myparameters: VecD::new(1),
            alldone: false,
            tolreached: false,
            my_multi_line_nb: 0,
            my_is_clear: false,
        }
    }

    /// OCCT Init (L1692-1709).
    pub fn init(
        &mut self,
        degreemin: usize,
        degreemax: usize,
        tolerance3d: f64,
        tolerance2d: f64,
        nb_iterations: i32,
        cutting: bool,
        parametrization: ApproxParamType,
    ) {
        self.mydegremin = degreemin;
        self.mydegremax = degreemax;
        self.mytol3d = tolerance3d;
        self.mytol2d = tolerance2d;
        self.par = parametrization;
        self.mycut = cutting;
        self.myitermax = nb_iterations;
    }

    pub fn set_constraints(&mut self, first_c: AppParConstraint, last_c: AppParConstraint) {
        self.myfirstc = first_c;
        self.mylastc = last_c;
    }

    pub fn is_all_approximated(&self) -> bool {
        self.alldone
    }
    pub fn is_tolerance_reached(&self) -> bool {
        self.tolreached
    }
    pub fn error(&self, index: usize) -> (f64, f64) {
        (self.tolers3d[index - 1], self.tolers2d[index - 1])
    }
    pub fn nb_multi_curves(&self) -> usize {
        self.mymulti_curves.len()
    }
    pub fn value(&self, index: usize) -> &MultiCurve {
        &self.mymulti_curves[index - 1]
    }
    pub fn value_mut(&mut self, index: usize) -> &mut MultiCurve {
        &mut self.mymulti_curves[index - 1]
    }

    /// OCCT Approx_ComputeLine::Parameters (L1249-1322) — chord-length /
    /// centripetal / iso-parametric parameters of the part points.
    fn parameters(&self, ml: &WLineAccess, first_p: usize, last_p: usize) -> VecD {
        approx_parameters(ml, first_p, last_p, self.par)
    }

    /// OCCT Approx_ComputeLine::FirstTangencyVector (L430-512).
    fn first_tangency_vector(&self, ml: &WLineAccess, index: usize) -> VecD {
        let nb_p3d = ml.nb_p3d();
        let nb_p2d = ml.nb_p2d();
        let dim = nb_p3d * 3 + nb_p2d * 2;
        let mut v = VecD::new(dim);
        if let Some((tab_v, tab_v2d)) = ml.tangency(index) {
            let mut j = 1usize;
            for i in 0..nb_p3d {
                v.set(j, tab_v[i].x);
                v.set(j + 1, tab_v[i].y);
                v.set(j + 2, tab_v[i].z);
                j += 3;
            }
            j = nb_p3d * 3 + 1;
            for i in 0..nb_p2d {
                v.set(j, tab_v2d[i].x);
                v.set(j + 1, tab_v2d[i].y);
                j += 2;
            }
        } else {
            // Parabola through index, index+1, index+2 (OCCT L481-511).
            let par = self.parameters(ml, index, index + 2);
            let ls = LeastSquare::new(
                ml,
                index,
                index + 2,
                AppParConstraint::PassPoint,
                AppParConstraint::PassPoint,
                3,
            );
            let mut ls = ls;
            ls.perform(&par);
            let c = ls.bezier_value();
            let (_, myv) = c.d1(1, 0.0);
            let mut j = 1usize;
            for i in 0..nb_p3d {
                v.set(j, myv.x);
                v.set(j + 1, myv.y);
                v.set(j + 2, myv.z);
                j += 3;
            }
            j = nb_p3d * 3 + 1;
            for i in 0..nb_p2d {
                let (_, myv2) = c.d1(nb_p3d + i + 1, 0.0);
                v.set(j, myv2.x);
                v.set(j + 1, myv2.y);
                j += 2;
            }
        }
        v
    }

    /// OCCT Approx_ComputeLine::LastTangencyVector (L514-595).
    fn last_tangency_vector(&self, ml: &WLineAccess, index: usize) -> VecD {
        let nb_p3d = ml.nb_p3d();
        let nb_p2d = ml.nb_p2d();
        let dim = nb_p3d * 3 + nb_p2d * 2;
        let mut v = VecD::new(dim);
        if let Some((tab_v, tab_v2d)) = ml.tangency(index) {
            let mut j = 1usize;
            for i in 0..nb_p3d {
                v.set(j, tab_v[i].x);
                v.set(j + 1, tab_v[i].y);
                v.set(j + 2, tab_v[i].z);
                j += 3;
            }
            j = nb_p3d * 3 + 1;
            for i in 0..nb_p2d {
                v.set(j, tab_v2d[i].x);
                v.set(j + 1, tab_v2d[i].y);
                j += 2;
            }
        } else {
            // Parabola through index-2, index-1, index (OCCT L564-593).
            let par = self.parameters(ml, index - 2, index);
            let ls = LeastSquare::new(
                ml,
                index - 2,
                index,
                AppParConstraint::PassPoint,
                AppParConstraint::PassPoint,
                3,
            );
            let mut ls = ls;
            ls.perform(&par);
            let c = ls.bezier_value();
            let mut j = 1usize;
            for i in 0..nb_p3d {
                let (_, myv) = c.d1(i + 1, 1.0);
                v.set(j, myv.x);
                v.set(j + 1, myv.y);
                v.set(j + 2, myv.z);
                j += 3;
            }
            j = nb_p3d * 3 + 1;
            for i in 0..nb_p2d {
                let (_, myv2) = c.d1(nb_p3d + i + 1, 1.0);
                v.set(j, myv2.x);
                v.set(j + 1, myv2.y);
                j += 2;
            }
        }
        v
    }

    /// OCCT Approx_ComputeLine::SearchFirstLambda (L597-655).
    fn search_first_lambda(&self, ml: &WLineAccess, the_param: &VecD, v: &VecD, index: usize) -> f64 {
        let nb_p3d = ml.nb_p3d();
        let nb_p2d = ml.nb_p2d();
        let p1 = ml.value_p3d(index);
        let p2 = ml.value_p3d(index + 1);
        let u1 = the_param.get(1);
        let u2 = the_param.get(2);
        let lambda;
        let s;
        if nb_p3d != 0 {
            let p1p2 = p2 - p1;
            let myv = DVec3::new(v.get(1), v.get(2), v.get(3));
            lambda = p1p2.length() / (myv.length() * (u2 - u1));
            s = if p1p2.dot(myv) > 0.0 { 1.0 } else { -1.0 };
        } else {
            let p12d = ml.value_p2d(index)[0];
            let p22d = ml.value_p2d(index + 1)[0];
            let p1p2 = p22d - p12d;
            let myv = DVec2::new(v.get(1), v.get(2));
            lambda = p1p2.length() / (myv.length() * (u2 - u1));
            s = if p1p2.dot(myv) > 0.0 { 1.0 } else { -1.0 };
        }
        let _ = nb_p2d;
        s * lambda
    }

    /// OCCT Approx_ComputeLine::SearchLastLambda (L657-715).
    fn search_last_lambda(&self, ml: &WLineAccess, the_param: &VecD, v: &VecD, index: usize) -> f64 {
        let nb_p3d = ml.nb_p3d();
        let nb_p2d = ml.nb_p2d();
        let p1 = ml.value_p3d(index - 1);
        let p2 = ml.value_p3d(index);
        let u1 = the_param.get(1);
        let u2 = the_param.get(2);
        let lambda;
        let s;
        if nb_p3d != 0 {
            let p1p2 = p2 - p1;
            let myv = DVec3::new(v.get(1), v.get(2), v.get(3));
            lambda = p1p2.length() / (myv.length() * (u2 - u1));
            s = if p1p2.dot(myv) > 0.0 { 1.0 } else { -1.0 };
        } else {
            let p12d = ml.value_p2d(index - 1)[0];
            let p22d = ml.value_p2d(index)[0];
            let p1p2 = p22d - p12d;
            let myv = DVec2::new(v.get(1), v.get(2));
            lambda = p1p2.length() / (myv.length() * (u2 - u1));
            s = if p1p2.dot(myv) > 0.0 { 1.0 } else { -1.0 };
        }
        let _ = nb_p2d;
        s * lambda
    }

    /// OCCT Approx_ComputeLine::Compute (L1324-1441).
    fn compute(
        &mut self,
        ml: &WLineAccess,
        fpt: usize,
        lpt: usize,
        para: &mut VecD,
        the_tol3d: &mut f64,
        the_tol2d: &mut f64,
        indbad: &mut usize,
    ) -> bool {
        *indbad = 0;
        let nbp = lpt - fpt + 1;
        let par_sav = para.clone();
        let mut m_degmax = self.mydegremax;
        if nbp < m_degmax + 5 && self.mycut {
            m_degmax = nbp - 5;
        }
        if m_degmax < self.mydegremin {
            m_degmax = self.mydegremin;
        }
        self.currenttol3d = f64::INFINITY;
        self.currenttol2d = f64::INFINITY;
        let deg_lo = (nbp - 1).min(self.mydegremin);
        for deg in deg_lo..=m_degmax {
            let grad = Gradient::new(
                ml,
                fpt,
                lpt,
                &self.myconstraints,
                para,
                deg,
                self.mytol3d,
                self.mytol2d,
                self.myitermax,
            );
            let mydone = grad.is_done();
            let my_scu = grad.value();
            if my_scu.nb_curves() == 0 {
                continue;
            }
            *the_tol3d = grad.max_error3d();
            *the_tol2d = grad.max_error2d();
            if std::env::var("RCAD_APX_DEBUG").is_ok() {
                eprintln!("  [APX] deg={} nbp={} mydone={} tol3d={:.3e} tol2d={:.3e} mytol3d={:.3e} mytol2d={:.3e} itermax={}", deg, nbp, mydone, *the_tol3d, *the_tol2d, self.mytol3d, self.mytol2d, self.myitermax);
            }
            // restau: restore parameters if not strictly increasing.
            let mut restau = false;
            let mut uu1 = para.get(1);
            for i in 2..=para.len() {
                let uu2 = para.get(i);
                if uu2 <= uu1 {
                    restau = true;
                    break;
                }
                uu1 = uu2;
            }
            if restau {
                for i in 1..=para.len() {
                    para.set(i, par_sav.get(i));
                }
            }
            if mydone {
                if *the_tol3d <= self.mytol3d && *the_tol2d <= self.mytol2d {
                    self.tolreached = true;
                    if !check_multicurve(&my_scu, ml, fpt, lpt, indbad) {
                        return false;
                    } else {
                        self.mymulti_curves.push(my_scu);
                        self.mypar.push(para.clone());
                        self.tolers3d.push(*the_tol3d);
                        self.tolers2d.push(*the_tol2d);
                        return true;
                    }
                }
            }
            if *the_tol3d <= self.currenttol3d && *the_tol2d <= self.currenttol2d {
                self.the_multi_curve = my_scu;
                self.currenttol3d = *the_tol3d;
                self.currenttol2d = *the_tol2d;
                self.myparameters = para.clone();
            }
        }
        false
    }

    /// OCCT Approx_ComputeLine::ComputeCurve (L1443-1690): the interpolation
    /// path used when the part has too few points.
    fn compute_curve(&mut self, ml: &WLineAccess, firstpt: usize, lastpt: usize) -> bool {
        let myfirstpt = firstpt;
        let mylastpt = lastpt;
        let nbp = lastpt - firstpt + 1;
        let nb_p3d = ml.nb_p3d();
        let nb_p2d = ml.nb_p2d();
        let para = self.parameters(ml, firstpt, lastpt);
        if nbp == 2 {
            // Linear interpolation poles (OCCT L1494-1642).
            let deg = self.mydegremin;
            let mut my_scu = MultiCurve::new(deg + 1, nb_p3d, nb_p2d);
            let p1 = ml.value_p3d(myfirstpt);
            let p2 = ml.value_p3d(mylastpt);
            let p2d1 = ml.value_p2d(myfirstpt);
            let p2d2 = ml.value_p2d(mylastpt);
            let mut mp1 = MultiPoint::new(nb_p3d, nb_p2d);
            for j in 0..nb_p3d {
                mp1.p3d[j] = p1;
            }
            for j in 0..nb_p2d {
                mp1.p2d[j] = p2d1[j];
            }
            let mut mp2 = MultiPoint::new(nb_p3d, nb_p2d);
            for j in 0..nb_p3d {
                mp2.p3d[j] = p2;
            }
            for j in 0..nb_p2d {
                mp2.p2d[j] = p2d2[j];
            }
            my_scu.set_value(1, mp1);
            my_scu.set_value(deg + 1, mp2);
            for i in 2..=deg {
                let t = (i - 1) as f64 / deg as f64;
                let mut mp = MultiPoint::new(nb_p3d, nb_p2d);
                for j in 0..nb_p3d {
                    mp.p3d[j] = p1 + (p2 - p1) * t;
                }
                for j in 0..nb_p2d {
                    mp.p2d[j] = p2d1[j] + (p2d2[j] - p2d1[j]) * t;
                }
                my_scu.set_value(i, mp);
            }
            self.tolreached = true;
            self.mymulti_curves.push(my_scu);
            self.mypar.push(para);
            self.tolers3d.push(1.0e-7);
            self.tolers2d.push(1.0e-10);
            return true;
        }
        // With the tangents (OCCT L1645-1689).
        let deg = nbp + 1;
        let mut my_scu = MultiCurve::new(deg + 1, nb_p3d, nb_p2d);
        let v1 = self.first_tangency_vector(ml, myfirstpt);
        let lambda1 = self.search_first_lambda(ml, &para, &v1, myfirstpt);
        let v2 = self.last_tangency_vector(ml, mylastpt);
        let lambda2 = self.search_last_lambda(ml, &para, &v2, mylastpt);
        let cons = AppParConstraint::TangencyPoint;
        let mut ls = LeastSquare::new(ml, myfirstpt, mylastpt, cons, cons, deg + 1);
        ls.vec1t = v1.v.clone();
        ls.vec2t = v2.v.clone();
        ls.first_constraint = AppParConstraint::TangencyPoint;
        ls.last_constraint = AppParConstraint::TangencyPoint;
        let l1 = lambda1 / deg as f64;
        let l2 = lambda2 / deg as f64;
        ls.lambda1 = l1;
        ls.lambda2 = l2;
        ls.tangency_fixed = true;
        let mut para2 = para.clone();
        ls.perform(&para2);
        let mydone = ls.is_done();
        my_scu = ls.bezier_value();
        if mydone {
            let (fv, the_tol3d, the_tol2d) = ls.error();
            let _ = fv;
            self.tolreached = true;
            self.mymulti_curves.push(my_scu);
            self.mypar.push(para2);
            self.tolers3d.push(the_tol3d);
            self.tolers2d.push(the_tol2d);
            return true;
        }
        false
    }

    /// OCCT Approx_ComputeLine::Perform (L832-1219) — the cutting loop.
    pub fn perform(&mut self, ml: &WLineAccess) {
        if !self.my_is_clear {
            self.mymulti_curves.clear();
            self.mypar.clear();
            self.tolers3d.clear();
            self.tolers2d.clear();
            self.my_multi_line_nb = 0;
        } else {
            self.my_is_clear = false;
        }
        let the_firstpt = ml.first_point();
        let the_lastpt = ml.last_point();
        let mut myfirstpt = the_firstpt;
        let mut mylastpt = the_lastpt;
        let mut finish = false;
        let mut begin = true;
        let mut ok = false;
        let mut go_up = false;
        let mut the_tol3d = 0.0;
        let mut the_tol2d = 0.0;
        let mut the_param = VecD::new(the_lastpt - the_firstpt + 1);
        self.myconstraints[0] = ConstraintCouple { index: myfirstpt as i32, constraint: self.myfirstc };
        self.myconstraints[1] = ConstraintCouple { index: mylastpt as i32, constraint: self.mylastc };
        while !finish {
            let mut oldlastpt = mylastpt;
            if !begin {
                if !go_up {
                    if ok {
                        myfirstpt = mylastpt;
                        mylastpt = the_lastpt;
                        if myfirstpt == the_lastpt {
                            finish = true;
                            self.alldone = true;
                            return;
                        }
                    } else {
                        let nbp = mylastpt - myfirstpt + 1;
                        let my_status = ml.what_status();
                        if my_status == ApproxStatus::NoPointsAdded && nbp <= self.mydegremax + 1 {
                            let interpol = self.compute_curve(ml, myfirstpt, mylastpt);
                            if interpol {
                                if mylastpt == the_lastpt {
                                    finish = true;
                                    self.alldone = true;
                                    return;
                                }
                            }
                        }
                        mylastpt = (myfirstpt + mylastpt) / 2;
                    }
                }
                go_up = false;
            }
            let nbp = mylastpt - myfirstpt + 1;
            let my_status = ml.what_status();
            if nbp <= self.mydegremax + 5 {
                go_up = false;
                ok = true;
                if my_status == ApproxStatus::PointsAdded {
                    // OCCT L974-1115: MakeMLBetween; on failure the
                    // parameterization switches to IsoParametric (L996) and the
                    // part is re-fitted (L997-1103).
                    go_up = true;
                    let an_other_line1 = ml.make_ml_between(myfirstpt, mylastpt, nbp - 1);
                    let nbpdsotherligne: i64 = match &an_other_line1 {
                        Some(l) => l.first_point() as i64 - l.last_point() as i64,
                        None => 0,
                    };
                    if nbpdsotherligne == 0 || self.my_multi_line_nb >= 3 {
                        // OCCT L985-1103: MakeML failed — fit with
                        // IsoParametric parameters.
                        if myfirstpt == mylastpt {
                            break;
                        }
                        self.myconstraints[0] = ConstraintCouple {
                            index: myfirstpt as i32,
                            constraint: self.myfirstc,
                        };
                        self.myconstraints[1] = ConstraintCouple {
                            index: mylastpt as i32,
                            constraint: self.mylastc,
                        };
                        let mut param = VecD::new(mylastpt - myfirstpt + 1);
                        let save_par = self.par;
                        self.par = ApproxParamType::IsoParametric;
                        let p = self.parameters(ml, myfirstpt, mylastpt);
                        for i in 1..=p.len() {
                            param.set(i, p.get(i));
                        }
                        self.the_multi_curve = MultiCurve::new(2, ml.nb_p3d(), ml.nb_p2d());
                        let mut an_other_line2: Option<WLineAccess> = None;
                        let mut is_other_line2_made = false;
                        let mut indbad = 0usize;
                        ok = self.compute(
                            ml,
                            myfirstpt,
                            mylastpt,
                            &mut param,
                            &mut the_tol3d,
                            &mut the_tol2d,
                            &mut indbad,
                        );
                        if indbad != 0 {
                            an_other_line2 = ml.make_ml_one_more_point(myfirstpt, mylastpt, indbad);
                            is_other_line2_made = an_other_line2.is_some();
                        }
                        if is_other_line2_made {
                            self.my_is_clear = true;
                            self.par = save_par;
                            if let Some(line2) = an_other_line2 {
                                self.perform(&line2);
                            }
                            ok = true;
                        }
                        if !ok {
                            // OCCT L1017-1051: ChordLength retry.
                            let tt3d = self.currenttol3d;
                            let tt2d = self.currenttol2d;
                            let save_parameters = self.myparameters.clone();
                            let save_multi_curve = self.the_multi_curve.clone();
                            if save_par != ApproxParamType::IsoParametric {
                                self.par = save_par;
                            } else {
                                self.par = ApproxParamType::ChordLength;
                            }
                            let p = self.parameters(ml, myfirstpt, mylastpt);
                            for i in 1..=p.len() {
                                param.set(i, p.get(i));
                            }
                            an_other_line2 = None;
                            is_other_line2_made = false;
                            indbad = 0;
                            ok = self.compute(
                                ml,
                                myfirstpt,
                                mylastpt,
                                &mut param,
                                &mut the_tol3d,
                                &mut the_tol2d,
                                &mut indbad,
                            );
                            if indbad != 0 {
                                an_other_line2 =
                                    ml.make_ml_one_more_point(myfirstpt, mylastpt, indbad);
                                is_other_line2_made = an_other_line2.is_some();
                            }
                            if is_other_line2_made {
                                self.my_is_clear = true;
                                if let Some(line2) = an_other_line2 {
                                    self.perform(&line2);
                                }
                                ok = true;
                            }
                            if !ok && tt3d <= self.currenttol3d && tt2d <= self.currenttol2d {
                                self.currenttol3d = tt3d;
                                self.currenttol2d = tt2d;
                                self.myparameters = save_parameters;
                                self.the_multi_curve = save_multi_curve;
                            }
                        }
                        self.par = save_par;
                        if myfirstpt == the_lastpt {
                            finish = true;
                            self.alldone = true;
                            return;
                        }
                        oldlastpt = mylastpt;
                        if !ok {
                            // OCCT L1062-1103: CheckMultiCurve + MakeMLOneMorePoint.
                            self.tolreached = false;
                            if self.the_multi_curve.nb_curves() == 0 {
                                self.mymulti_curves.clear();
                                return;
                            }
                            let mut an_other_line3: Option<WLineAccess> = None;
                            let mut indbad2 = 0usize;
                            if !check_multicurve(
                                &self.the_multi_curve,
                                ml,
                                myfirstpt,
                                mylastpt,
                                &mut indbad2,
                            ) {
                                an_other_line3 =
                                    ml.make_ml_one_more_point(myfirstpt, mylastpt, indbad2);
                            }
                            if let Some(line3) = an_other_line3 {
                                self.my_is_clear = true;
                                self.perform(&line3);
                                myfirstpt = mylastpt;
                                mylastpt = the_lastpt;
                            } else {
                                self.mymulti_curves.push(self.the_multi_curve.clone());
                                self.tolers3d.push(self.currenttol3d);
                                self.tolers2d.push(self.currenttol2d);
                                let mylen = oldlastpt - myfirstpt + 1;
                                let my_par_len = self.myparameters.len();
                                let a_len = my_par_len.max(mylen);
                                let mut the_par = VecD::new(a_len);
                                for i in 0..a_len {
                                    the_par.set(
                                        i + 1,
                                        if i < self.myparameters.len() {
                                            self.myparameters.get(i + 1)
                                        } else {
                                            0.0
                                        },
                                    );
                                }
                                self.mypar.push(the_par);
                                myfirstpt = oldlastpt;
                                mylastpt = the_lastpt;
                            }
                        }
                        // OCCT L1119-1120: advance to the next part.
                        myfirstpt = oldlastpt;
                        mylastpt = the_lastpt;
                    } else {
                        // OCCT L1108-1115: MakeML succeeded — recurse on the
                        // densified line.
                        self.my_is_clear = true;
                        self.my_multi_line_nb += 1;
                        if let Some(line1) = an_other_line1 {
                            self.perform(&line1);
                        }
                        myfirstpt = mylastpt;
                        mylastpt = the_lastpt;
                    }
                }
                // OCCT L1118-1147: NoPointsAdded with a small part — keep
                // the best approximation obtained so far.  This runs
                // BEFORE the Compute of this interval; GoUp then skips it
                // (so currenttol3d still holds the previous interval's
                // best effort — a finite value, never the Compute's
                // freshly-reset RealLast/inf).
                if my_status == ApproxStatus::NoPointsAdded && !begin {
                    go_up = true;
                    self.tolreached = false;
                    if self.the_multi_curve.nb_curves() == 0 {
                        self.mymulti_curves.clear();
                        return;
                    }
                    self.mymulti_curves.push(self.the_multi_curve.clone());
                    self.tolers3d.push(self.currenttol3d);
                    self.tolers2d.push(self.currenttol2d);
                    let mylen = oldlastpt - myfirstpt + 1;
                    let my_par_len = self.myparameters.len();
                    let a_len = my_par_len.max(mylen);
                    let mut the_par = VecD::new(a_len);
                    for i in 0..a_len {
                        the_par.set(i + 1, if i < self.myparameters.len() { self.myparameters.get(i + 1) } else { 0.0 });
                    }
                    self.mypar.push(the_par);
                    myfirstpt = oldlastpt;
                    mylastpt = the_lastpt;
                } else if my_status == ApproxStatus::NoApproximation {
                    // OCCT L1149-1157: no approximation is done between
                    // myfirstpt and mylastpt.
                    go_up = true;
                    myfirstpt = mylastpt;
                    mylastpt = the_lastpt;
                }
            }
            if myfirstpt == the_lastpt {
                finish = true;
                self.alldone = true;
                return;
            }
            if !go_up {
                if myfirstpt == mylastpt {
                    break;
                }
                self.myconstraints[0] = ConstraintCouple { index: myfirstpt as i32, constraint: self.myfirstc };
                self.myconstraints[1] = ConstraintCouple { index: mylastpt as i32, constraint: self.mylastc };
                let mut param = VecD::new(mylastpt - myfirstpt + 1);
                if begin {
                    if self.myfirst_param.is_none() {
                        let p = self.parameters(ml, myfirstpt, mylastpt);
                        for i in 1..=p.len() {
                            param.set(i, p.get(i));
                        }
                    } else {
                        let fp = self.myfirst_param.as_ref().unwrap();
                        for i in 1..=fp.len() {
                            param.set(i, fp.get(i));
                        }
                        self.myfirst_param = None;
                    }
                    for i in 1..=param.len() {
                        the_param.set(i, param.get(i));
                    }
                    begin = false;
                } else {
                    let pfirst = the_param.get(myfirstpt - the_firstpt + 1);
                    let plast = the_param.get(mylastpt - the_firstpt + 1);
                    for i in myfirstpt..=mylastpt {
                        let v = (the_param.get(i - the_firstpt + 1) - pfirst) / (plast - pfirst);
                        param.set(i - myfirstpt + 1, v);
                    }
                }
                self.the_multi_curve = MultiCurve::new(2, ml.nb_p3d(), ml.nb_p2d());
                let mut indbad = 0usize;
                ok = self.compute(ml, myfirstpt, mylastpt, &mut param, &mut the_tol3d, &mut the_tol2d, &mut indbad);
                if myfirstpt == the_lastpt {
                    finish = true;
                    self.alldone = true;
                    return;
                }
            }
        }
    }
}

// ============================================================================
// AppParCurves_LeastSquare — semantic equivalent of
// AppParCurves_LeastSquare.gxx (L138-734): constrained least-squares fit of
// a Bezier multi-curve.  The OCCT dense solves (math_Householder, DACTCL
// banded LDLT) are replaced by the normal equations with Gaussian
// elimination; the equation structure and the constraint handling follow the
// OCCT source 1:1.
// ============================================================================

/// OCCT AppParCurves_LeastSquare::TheFirstPoint (L1465-1472).
fn the_first_point(cons: AppParConstraint, first_point: usize) -> usize {
    if cons == AppParConstraint::NoConstraint {
        first_point
    } else {
        first_point + 1
    }
}

/// OCCT AppParCurves_LeastSquare::TheLastPoint (L1474-1481).
fn the_last_point(cons: AppParConstraint, last_point: usize) -> usize {
    if cons == AppParConstraint::NoConstraint {
        last_point
    } else {
        last_point - 1
    }
}

/// Number of coordinate columns per multipoint (nbP*3 + nbP2d*2).
fn nb_b_columns(nb_p: usize, nb_p2d: usize) -> usize {
    nb_p * 3 + nb_p2d * 2
}

pub struct LeastSquare {
    // OCCT members.
    pub myfirstp: usize,
    pub mylastp: usize,
    pub first_p: usize,
    pub last_p: usize,
    pub nbpoles: usize,
    pub deg: usize,
    pub nb_p: usize,
    pub nb_p2d: usize,
    pub n_bcols: usize,
    pub first_constraint: AppParConstraint,
    pub last_constraint: AppParConstraint,
    pub mypoints: Vec<Vec<f64>>, // [data point (myfirstp..mylastp)][bcol]
    pub vec1t: Vec<f64>,
    pub vec2t: Vec<f64>,
    pub a: Vec<Vec<f64>>, // [data point][pole] Bernstein basis
    pub da: Vec<Vec<f64>>, // derivative of the basis
    pub mypoles: Vec<Vec<f64>>, // [pole][bcol]
    pub done: bool,
    /// OCCT Perform(Parameters, V1t, V2t, l1, l2): preset tangency scalars
    /// (the ComputeCurve interpolation path); when set, the tangency poles
    /// are fixed instead of being free unknowns.
    pub tangency_fixed: bool,
    pub lambda1: f64,
    pub lambda2: f64,
    // OCCT LeastSquare members (Init L274-508).
    pub resinit: usize,
    pub resfin: usize,
    pub isready: bool,
    pub iscalculated: bool,
    pub na: usize,
    pub nlignes: usize,
    pub ninc: usize,
    /// OCCT B2: the reduced right-hand side rows (FirstP..LastP).
    pub b2: Vec<Vec<f64>>,
}

impl LeastSquare {
    /// OCCT LeastSquare(SSP, FirstPoint, LastPoint, FirstCons, LastCons,
    /// NbPol) constructor + Init (L138-167, L274-508).
    pub fn new(
        ml: &WLineAccess,
        first_point: usize,
        last_point: usize,
        mut first_cons: AppParConstraint,
        mut last_cons: AppParConstraint,
        nbpoles: usize,
    ) -> Self {
        let nb_p = ml.nb_p3d();
        let nb_p2d = ml.nb_p2d();
        let n_bcols = nb_b_columns(nb_p, nb_p2d);
        let mut ls = LeastSquare {
            myfirstp: first_point,
            mylastp: last_point,
            first_p: 0,
            last_p: 0,
            nbpoles,
            deg: nbpoles - 1,
            nb_p,
            nb_p2d,
            n_bcols,
            first_constraint: first_cons,
            last_constraint: last_cons,
            mypoints: vec![vec![0.0; n_bcols]; last_point - first_point + 1],
            vec1t: vec![0.0; n_bcols],
            vec2t: vec![0.0; n_bcols],
            a: vec![vec![0.0; nbpoles]; last_point - first_point + 1],
            da: vec![vec![0.0; nbpoles]; last_point - first_point + 1],
            mypoles: vec![vec![0.0; n_bcols]; nbpoles],
            done: false,
            tangency_fixed: false,
            lambda1: 0.0,
            lambda2: 0.0,
            resinit: 0,
            resfin: 0,
            isready: false,
            iscalculated: false,
            na: 0,
            nlignes: 0,
            ninc: 0,
            b2: vec![vec![0.0; n_bcols]; last_point.saturating_sub(first_point)],
        };
        ls.first_p = the_first_point(first_cons, first_point);
        ls.last_p = the_last_point(last_cons, last_point);
        // OCCT Affect (L1064-1205): tangent vectors; with no SvSurfaces the
        // TangencyPoint constraint degrades to PassPoint.
        ls.affect(first_point, &mut first_cons, ml);
        ls.affect(last_point, &mut last_cons, ml);
        ls.first_constraint = first_cons;
        ls.last_constraint = last_cons;
        // mypoints: the translated points.
        for j in 0..(last_point - first_point + 1) {
            let idx = first_point + j;
            let p3 = ml.value_p3d(idx);
            let p2 = ml.value_p2d(idx);
            let mut i2 = 0;
            for i in 0..nb_p {
                ls.mypoints[j][i2] = p3.x;
                ls.mypoints[j][i2 + 1] = p3.y;
                ls.mypoints[j][i2 + 2] = p3.z;
                i2 += 3;
            }
            for i in 0..nb_p2d {
                ls.mypoints[j][i2] = p2[i].x;
                ls.mypoints[j][i2 + 1] = p2[i].y;
                i2 += 2;
            }
        }
        // Fixed end poles.
        if first_cons != AppParConstraint::NoConstraint {
            for i in 0..n_bcols {
                ls.mypoles[0][i] = ls.mypoints[0][i];
            }
        }
        if last_cons != AppParConstraint::NoConstraint {
            for i in 0..n_bcols {
                ls.mypoles[nbpoles - 1][i] = ls.mypoints[last_point - first_point][i];
            }
        }
        // OCCT Init L396-508: free pole range, NA/Nlignes/Ninc.
        ls.resinit = match ls.first_constraint {
            AppParConstraint::NoConstraint => 1,
            AppParConstraint::PassPoint => 2,
            AppParConstraint::TangencyPoint => 3,
            _ => 4,
        };
        ls.resfin = match ls.last_constraint {
            AppParConstraint::NoConstraint => nbpoles,
            AppParConstraint::PassPoint => nbpoles - 1,
            AppParConstraint::TangencyPoint => nbpoles - 2,
            _ => nbpoles - 3,
        };
        let nincx = if ls.resfin >= ls.resinit { ls.resfin - ls.resinit + 1 } else { 0 };
        if nincx < 1 {
            ls.isready = false;
        } else {
            ls.isready = true;
            let neq = ls.last_p - ls.first_p + 1;
            ls.na = 3 * nb_p + 2 * nb_p2d;
            ls.nlignes = ls.na * neq;
            ls.ninc = ls.na * nincx;
            if ls.first_constraint >= AppParConstraint::TangencyPoint {
                ls.ninc += 1;
            }
            if ls.last_constraint >= AppParConstraint::TangencyPoint {
                ls.ninc += 1;
            }
            ls.b2 = vec![vec![0.0; n_bcols]; neq];
        }
        ls
    }

    /// OCCT Affect (L1064-1205): fill the tangent (Vt) vectors at `index`;
    /// when the line has no tangency the constraint degrades to PassPoint.
    fn affect(&mut self, index: usize, cons: &mut AppParConstraint, ml: &WLineAccess) {
        if *cons >= AppParConstraint::TangencyPoint {
            let vt = if index == self.myfirstp {
                &mut self.vec1t
            } else {
                &mut self.vec2t
            };
            if let Some((mut v3, mut v2d)) = ml.tangency(index) {
                // OCCT CheckTangents (AppParCurves_LeastSquare.gxx L64-106):
                // the tangent direction must agree with the direction between
                // the adjacent points; otherwise reverse ALL tangents.
                let (pt1, pt2) = if index < ml.last_point() {
                    (ml.value_p3d(index), ml.value_p3d(index + 1))
                } else {
                    (ml.value_p3d(index - 1), ml.value_p3d(index))
                };
                let mut is_to_change_dir = false;
                for i in 0..self.nb_p {
                    let a_v1 = pt2 - pt1;
                    let a_v2 = v3[i];
                    if a_v1.dot(a_v2) < 0.0 {
                        is_to_change_dir = true;
                        break;
                    }
                }
                if is_to_change_dir {
                    for v in v3.iter_mut() {
                        *v = -*v;
                    }
                    for v in v2d.iter_mut() {
                        *v = -*v;
                    }
                }
                let mut i2 = 0;
                for i in 0..self.nb_p {
                    vt[i2] = v3[i].x;
                    vt[i2 + 1] = v3[i].y;
                    vt[i2 + 2] = v3[i].z;
                    i2 += 3;
                }
                for i in 0..self.nb_p2d {
                    vt[i2] = v2d[i].x;
                    vt[i2 + 1] = v2d[i].y;
                    i2 += 2;
                }
            } else {
                *cons = AppParConstraint::PassPoint;
            }
        }
    }

    /// OCCT ComputeFunction (L1483-1493): the Bernstein basis and its
    /// derivative at the given parameters.
    pub fn compute_function(&mut self, parameters: &VecD) {
        for j in 0..(self.mylastp - self.myfirstp + 1) {
            let u = parameters.get(j + 1);
            for p in 0..self.nbpoles {
                let (b, db) = bernstein(self.deg, p, u);
                self.a[j][p] = b;
                self.da[j][p] = db;
            }
        }
    }

    /// OCCT Perform(Parameters) (L510-734): solve the constrained least
    /// squares.  The NoConstraint case uses math_Householder (L535-544); the
    /// PassPoint cases build the reduced system B2 and solve it with
    /// SearchIndex + MakeTAA + DACTCL (L557-620).  The tangency cases
    /// (L622-734) are not reachable for the quadric WLine: the TangencyPoint
    /// constraint degrades to PassPoint in Affect.
    pub fn perform(&mut self, parameters: &VecD) {
        self.done = false;
        if !self.isready {
            return;
        }
        self.iscalculated = false;
        self.compute_function(parameters);
        let n_bcols = self.n_bcols;
        let nbpoles = self.nbpoles;
        if self.first_constraint != AppParConstraint::TangencyPoint
            && self.last_constraint != AppParConstraint::TangencyPoint
        {
            if self.first_constraint == AppParConstraint::NoConstraint {
                if self.last_constraint == AppParConstraint::NoConstraint {
                    // OCCT L535-544: math_Householder HouResol(A, mypoints).
                    let a_mat = rcad_kernel::math::MatD { m: self.a.clone() };
                    let b_mat = rcad_kernel::math::MatD { m: self.mypoints.clone() };
                    match rcad_kernel::math::lin::householder_ls(&a_mat, &b_mat, 1.0e-20) {
                        Some(sol) => {
                            self.mypoles = sol.m;
                            self.done = true;
                            return;
                        }
                        None => {
                            self.done = false;
                            return;
                        }
                    }
                } else {
                    // OCCT L547-554: B2 = mypoints - A(:,nbpoles)*poleN.
                    for j in self.first_p..=self.last_p {
                        let dj = j - self.myfirstp;
                        for i in 0..n_bcols {
                            self.b2[j - self.first_p][i] = self.mypoints[dj][i]
                                - self.a[dj][nbpoles - 1] * self.mypoles[nbpoles - 1][i];
                        }
                    }
                }
            } else if self.first_constraint == AppParConstraint::PassPoint {
                if self.last_constraint == AppParConstraint::NoConstraint {
                    // OCCT L561-568: B2 = mypoints - A(:,1)*pole1.
                    for j in self.first_p..=self.last_p {
                        let dj = j - self.myfirstp;
                        for i in 0..n_bcols {
                            self.b2[j - self.first_p][i] =
                                self.mypoints[dj][i] - self.a[dj][0] * self.mypoles[0][i];
                        }
                    }
                } else if self.last_constraint == AppParConstraint::PassPoint {
                    // OCCT L572-580: B2 = mypoints - A(:,1)*pole1 - A(:,nbpoles)*poleN.
                    for j in self.first_p..=self.last_p {
                        let dj = j - self.myfirstp;
                        for i in 0..n_bcols {
                            self.b2[j - self.first_p][i] = self.mypoints[dj][i]
                                - self.a[dj][0] * self.mypoles[0][i]
                                - self.a[dj][nbpoles - 1] * self.mypoles[nbpoles - 1][i];
                        }
                    }
                }
            }
            // OCCT L586-620: resolution.
            let nincx = if self.resfin >= self.resinit {
                self.resfin - self.resinit + 1
            } else {
                0
            };
            if nincx < 1 {
                self.done = true;
                return;
            }
            let mut index = rcad_kernel::math::IntVec::new(nincx);
            self.search_index(&mut index);
            let mut mytab = rcad_kernel::math::MatD::new(nincx, n_bcols);
            let the_aa_len = index.get(nincx) as usize;
            let mut the_aa = rcad_kernel::math::VecD::new(the_aa_len);
            self.make_taa(&mut the_aa, &mut mytab);
            if rcad_kernel::math::lin::dactcl_decompose(&mut the_aa, &index.v, 1.0e-20) != 0 {
                self.done = false;
                return;
            }
            for j in 1..=n_bcols {
                let mut my_tabb = rcad_kernel::math::VecD::new(nincx);
                for (kk2, i) in (self.resinit..=self.resfin).enumerate() {
                    my_tabb.set(kk2 + 1, mytab.get(i - self.resinit + 1, j));
                }
                if rcad_kernel::math::lin::dactcl_solve(&the_aa, &mut my_tabb, &index.v, 1.0e-20) != 0 {
                    self.done = false;
                    return;
                }
                let mut i2 = 0usize;
                for k in self.resinit..=self.resfin {
                    self.mypoles[k - 1][j - 1] = my_tabb.get(i2 + 1);
                    i2 += 1;
                }
            }
            self.done = true;
        } else {
            // OCCT L622-734 (tangency constraints): the lambda parameters are
            // part of the DACTCL system (lambda1/lambda2 from the solved
            // myTAB) and the tangency poles are mypoints +- lambda*Vt.
            self.perform_tangency_full();
        }
    }

    /// OCCT Perform(Parameters) tangency branch (AppParCurves_LeastSquare.gxx
    /// L622-734): the tangency lambdas are unknowns of the DACTCL system
    /// (lambda1/lambda2 from the solved myTAB, L669-677) and the tangency
    /// poles are mypoints +- lambda*Vt (L695-706, L721-730).
    fn perform_tangency_full(&mut self) {
        let n_bcols = self.n_bcols;
        let nbpoles = self.nbpoles;
        let nincx = if self.resfin >= self.resinit {
            self.resfin - self.resinit + 1
        } else {
            0
        };
        let nincx2 = 2 * nincx;
        let ninc1 = if self.first_constraint >= AppParConstraint::TangencyPoint
            && self.last_constraint >= AppParConstraint::TangencyPoint
        {
            self.ninc - 1
        } else {
            self.ninc
        };
        let mut internal_index = rcad_kernel::math::IntVec::new(nincx);
        self.search_index(&mut internal_index);
        // The pivot index (L629-656).
        let mut index = rcad_kernel::math::IntVec::new(self.ninc);
        let mut l = 1usize;
        if self.resinit <= self.resfin {
            for j in 0..self.na {
                let deport = j * internal_index.get(nincx) as usize;
                for i in 1..=nincx {
                    index.set(l, internal_index.get(i) + deport as i32);
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
                index.set(ninc1, index.get(ninc1 - 1) + ninc1 as i32);
            }
        }
        if self.first_constraint >= AppParConstraint::TangencyPoint
            || self.last_constraint >= AppParConstraint::TangencyPoint
        {
            index.set(self.ninc, index.get(self.ninc - 1) + self.ninc as i32);
        }
        let mut the_a = rcad_kernel::math::VecD::new(index.get(self.ninc) as usize);
        let mut my_tab = rcad_kernel::math::VecD::new(self.ninc);
        self.make_taa_tangency(&mut the_a, &mut my_tab);
        let dec_err = rcad_kernel::math::lin::dactcl_decompose(&mut the_a, &index.v, 1.0e-20);
        if std::env::var("RCAD_TG_DEBUG").is_ok() {
            eprintln!("[TGD] dec_err={} index={:?} the_a={:?} my_tab={:?}", dec_err, index.v, &the_a.v[0..the_a.v.len().min(30)], &my_tab.v[0..my_tab.v.len().min(30)]);
        }
        if dec_err != 0 {
            self.done = false;
            return;
        }
        if rcad_kernel::math::lin::dactcl_solve(&the_a, &mut my_tab, &index.v, 1.0e-20) != 0 {
            self.done = false;
            return;
        }
        if std::env::var("RCAD_TG_DEBUG").is_ok() {
            eprintln!("[TGDS] solved={:?}", &my_tab.v[0..my_tab.v.len().min(30)]);
        }
        self.done = true;
        // lambda1/lambda2 (L669-677).
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
        if std::env::var("RCAD_TG_DEBUG").is_ok() {
            eprintln!("[TGS] ninc={} nincx={} ninc1={} l1={:.6e} l2={:.6e} v1t={:?} v2t={:?}", self.ninc, nincx, ninc1, self.lambda1, self.lambda2, &self.vec1t[0..3.min(self.vec1t.len())], &self.vec2t[0..3.min(self.vec2t.len())]);
        }
        // The mypoles fill (L681-733).
        let last_rel = self.mylastp - self.myfirstp;
        let mut k = 0usize;
        let mut i2 = 1usize;
        for _ci in 0..self.nb_p {
            for j in self.resinit..=self.resfin {
                let p = j - 1;
                self.mypoles[p][k] = my_tab.get(i2);
                self.mypoles[p][k + 1] = my_tab.get(i2 + nincx);
                self.mypoles[p][k + 2] = my_tab.get(i2 + nincx2);
                i2 += 1;
            }
            if self.first_constraint >= AppParConstraint::TangencyPoint {
                self.mypoles[1][k] = self.mypoints[0][k] + self.lambda1 * self.vec1t[k];
                self.mypoles[1][k + 1] = self.mypoints[0][k + 1] + self.lambda1 * self.vec1t[k + 1];
                self.mypoles[1][k + 2] = self.mypoints[0][k + 2] + self.lambda1 * self.vec1t[k + 2];
            }
            if self.last_constraint >= AppParConstraint::TangencyPoint {
                self.mypoles[nbpoles - 2][k] =
                    self.mypoints[last_rel][k] - self.lambda2 * self.vec2t[k];
                self.mypoles[nbpoles - 2][k + 1] =
                    self.mypoints[last_rel][k + 1] - self.lambda2 * self.vec2t[k + 1];
                self.mypoles[nbpoles - 2][k + 2] =
                    self.mypoints[last_rel][k + 2] - self.lambda2 * self.vec2t[k + 2];
            }
            k += 3;
            i2 += nincx2;
        }
        for _ci in 0..self.nb_p2d {
            for j in self.resinit..=self.resfin {
                let p = j - 1;
                self.mypoles[p][k] = my_tab.get(i2);
                self.mypoles[p][k + 1] = my_tab.get(i2 + nincx);
                i2 += 1;
            }
            if self.first_constraint >= AppParConstraint::TangencyPoint {
                self.mypoles[1][k] = self.mypoints[0][k] + self.lambda1 * self.vec1t[k];
                self.mypoles[1][k + 1] = self.mypoints[0][k + 1] + self.lambda1 * self.vec1t[k + 1];
            }
            if self.last_constraint >= AppParConstraint::TangencyPoint {
                self.mypoles[nbpoles - 2][k] =
                    self.mypoints[last_rel][k] - self.lambda2 * self.vec2t[k];
                self.mypoles[nbpoles - 2][k + 1] =
                    self.mypoints[last_rel][k + 1] - self.lambda2 * self.vec2t[k + 1];
            }
            k += 2;
            i2 += nincx;
        }
        if std::env::var("RCAD_TG_DEBUG").is_ok() {
            let p0 = &self.mypoles[0];
            let p1 = &self.mypoles[1];
            let pm = &self.mypoles[nbpoles / 2];
            let pn2 = &self.mypoles[nbpoles - 2];
            let pn = &self.mypoles[nbpoles - 1];
            eprintln!("[TGSP] P0=({:.6},{:.6},{:.6}) P1=({:.6},{:.6},{:.6}) Pm=({:.6},{:.6},{:.6}) Pn2=({:.6},{:.6},{:.6}) Pn=({:.6},{:.6},{:.6})", p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], pm[0], pm[1], pm[2], pn2[0], pn2[1], pn2[2], pn[0], pn[1], pn[2]);
        }
    }

    /// OCCT LeastSquare::MakeTAA(TheA, myTAB) (AppParCurves_LeastSquare.gxx
    /// L1510-1703): the full normal equations including the tangency (lambda)
    /// rows.  The interior block is the per-coordinate copy of the reduced
    /// triangle (MakeTAA(AA) L1763-1813) and the tangency rows are appended.
    fn make_taa_tangency(&self, the_a: &mut rcad_kernel::math::VecD, my_tab: &mut rcad_kernel::math::VecD) {
        let nincx = if self.resfin >= self.resinit {
            self.resfin - self.resinit + 1
        } else {
            0
        };
        let neq = self.last_p - self.first_p + 1;
        let ninc1 = if self.first_constraint >= AppParConstraint::TangencyPoint
            && self.last_constraint >= AppParConstraint::TangencyPoint
        {
            self.ninc - 1
        } else {
            self.ninc
        };
        let na1 = self.na - 1;
        let nlignes = self.nlignes;
        let mut my_b = vec![0.0f64; nlignes];
        let mut my_v1 = vec![0.0f64; nlignes];
        let mut my_v2 = vec![0.0f64; nlignes];
        let mut the_v1 = vec![0.0f64; self.ninc];
        let mut the_v2 = vec![0.0f64; self.ninc];
        let mut taf1 = 0.0;
        let mut taf2 = 0.0;
        let mut taf3 = 0.0;
        let mut tab1 = 0.0;
        let mut tab2 = 0.0;
        let last_rel = self.mylastp - self.myfirstp;
        if std::env::var("RCAD_TG_DEBUG").is_ok() {
            let di0 = self.first_p - self.myfirstp;
            eprintln!("[B] myfirstp={} mylastp={} first_p={} last_p={} di0={} last_rel={} mypoints0={:?} mypointsN={:?}", self.myfirstp, self.mylastp, self.first_p, self.last_p, di0, last_rel, &self.mypoints[di0][0..5], &self.mypoints[last_rel][0..5]);
            eprintln!("[B] a0={:?} aN={:?} vec1t={:?} vec2t={:?}", &self.a[di0][0..3], &self.a[last_rel][0..3], &self.vec1t[0..3], &self.vec2t[0..3]);
        }
        // myB/myV1/myV2 (L1540-1600).
        for i in self.first_p..=self.last_p {
            let di = i - self.myfirstp;
            let ai2 = self.a[di][1];
            let aid = self.a[di][self.nbpoles - 2];
            if std::env::var("RCAD_TG_DEBUG").is_ok() && i == self.first_p {
                let xx0 = if self.first_constraint >= AppParConstraint::PassPoint { self.a[di][0] } else { 0.0 };
                let yy0 = if self.last_constraint >= AppParConstraint::PassPoint { self.a[di][self.nbpoles - 1] } else { 0.0 };
                eprintln!("[BB] i={} xx={:.6} yy={:.6} bxyz=({:.6},{:.6},{:.6}) mypoints=({:.6},{:.6},{:.6}) Pn=({:.6},{:.6},{:.6})", i, xx0 + ai2, yy0 + aid, self.mypoints[di][0] - (xx0 + ai2) * self.mypoints[0][0] - (yy0 + aid) * self.mypoints[last_rel][0], self.mypoints[di][1] - (xx0 + ai2) * self.mypoints[0][1] - (yy0 + aid) * self.mypoints[last_rel][1], self.mypoints[di][2] - (xx0 + ai2) * self.mypoints[0][2] - (yy0 + aid) * self.mypoints[last_rel][2], self.mypoints[di][0], self.mypoints[di][1], self.mypoints[di][2], self.mypoints[last_rel][0], self.mypoints[last_rel][1], self.mypoints[last_rel][2]);
            }
            let mut xx = 0.0;
            let mut yy = 0.0;
            if self.first_constraint >= AppParConstraint::PassPoint {
                xx = self.a[di][0];
            }
            if self.first_constraint >= AppParConstraint::TangencyPoint {
                xx += ai2;
            }
            if self.last_constraint >= AppParConstraint::PassPoint {
                yy = self.a[di][self.nbpoles - 1];
            }
            if self.last_constraint >= AppParConstraint::TangencyPoint {
                yy += aid;
            }
            let mut i2 = 0usize;
            let mut nrow = 0usize;
            let ix0 = i - self.first_p;
            for _ci in 0..self.nb_p {
                let ix = ix0 + nrow;
                let iy = ix + neq;
                let iz = iy + neq;
                if self.first_constraint >= AppParConstraint::TangencyPoint {
                    my_v1[ix] = ai2 * self.vec1t[i2];
                    my_v1[iy] = ai2 * self.vec1t[i2 + 1];
                    my_v1[iz] = ai2 * self.vec1t[i2 + 2];
                }
                if self.last_constraint >= AppParConstraint::TangencyPoint {
                    my_v2[ix] = -aid * self.vec2t[i2];
                    my_v2[iy] = -aid * self.vec2t[i2 + 1];
                    my_v2[iz] = -aid * self.vec2t[i2 + 2];
                }
                my_b[ix] = self.mypoints[di][i2] - xx * self.mypoints[0][i2]
                    - yy * self.mypoints[last_rel][i2];
                my_b[iy] = self.mypoints[di][i2 + 1] - xx * self.mypoints[0][i2 + 1]
                    - yy * self.mypoints[last_rel][i2 + 1];
                my_b[iz] = self.mypoints[di][i2 + 2] - xx * self.mypoints[0][i2 + 2]
                    - yy * self.mypoints[last_rel][i2 + 2];
                i2 += 3;
                nrow += 3 * neq;
            }
            for _ci in 0..self.nb_p2d {
                let ix = ix0 + nrow;
                let iy = ix + neq;
                if self.first_constraint >= AppParConstraint::TangencyPoint {
                    my_v1[ix] = ai2 * self.vec1t[i2];
                    my_v1[iy] = ai2 * self.vec1t[i2 + 1];
                }
                if self.last_constraint >= AppParConstraint::TangencyPoint {
                    my_v2[ix] = -aid * self.vec2t[i2];
                    my_v2[iy] = -aid * self.vec2t[i2 + 1];
                }
                my_b[ix] = self.mypoints[di][i2] - xx * self.mypoints[0][i2]
                    - yy * self.mypoints[last_rel][i2];
                my_b[iy] = self.mypoints[di][i2 + 1] - xx * self.mypoints[0][i2 + 1]
                    - yy * self.mypoints[last_rel][i2 + 1];
                nrow += 2 * neq;
                i2 += 2;
            }
        }
        // The normal equations (L1605-1647).
        for k in self.first_p..=self.last_p {
            let dk = k - self.myfirstp;
            let jinit = self.resinit;
            let jfin = self.resfin.min(self.nbpoles);
            let k1 = k - self.first_p;
            for i in 0..=na1 {
                let nb = i * neq + k1;
                let mut v1 = 0.0;
                let mut v2 = 0.0;
                if self.first_constraint >= AppParConstraint::TangencyPoint {
                    v1 = my_v1[nb];
                }
                if self.last_constraint >= AppParConstraint::TangencyPoint {
                    v2 = my_v2[nb];
                }
                let b = my_b[nb];
                let inc = i as i64 * nincx as i64 + 1 - self.resinit as i64;
                for j in jinit..=jfin {
                    let akj = self.a[dk][j - 1];
                    let u = (j as i64 + inc - 1) as usize;
                    if self.first_constraint >= AppParConstraint::TangencyPoint {
                        the_v1[u] += akj * v1;
                    }
                    if self.last_constraint >= AppParConstraint::TangencyPoint {
                        the_v2[u] += akj * v2;
                    }
                    my_tab.v[u] += akj * b;
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
        // The lambda diagonal (L1649-1662).
        if self.first_constraint >= AppParConstraint::TangencyPoint {
            the_v1[ninc1 - 1] = taf1;
            my_tab.set(ninc1, tab1);
        }
        if self.last_constraint >= AppParConstraint::TangencyPoint {
            the_v2[self.ninc - 1] = taf2;
            my_tab.set(self.ninc, tab2);
        }
        if self.first_constraint >= AppParConstraint::TangencyPoint
            && self.last_constraint >= AppParConstraint::TangencyPoint
        {
            the_v2[ninc1 - 1] = taf3;
        }
        // The interior block (L1664-1680): the per-coordinate copies of the
        // reduced triangle.
        if self.resinit <= self.resfin {
            let mut index = rcad_kernel::math::IntVec::new(nincx);
            self.search_index(&mut index);
            let mut aa = rcad_kernel::math::VecD::new(index.get(nincx) as usize);
            let mut my_tab2 = rcad_kernel::math::MatD::new(nincx, self.n_bcols);
            self.make_taa(&mut aa, &mut my_tab2);
            let mut kk = 0usize;
            for _k in 0..self.na {
                for i in 0..aa.len() {
                    the_a.set(kk + 1, aa.v[i]);
                    kk += 1;
                }
            }
        }
        // The tangency rows (L1684-1702).
        let length = the_a.len();
        if self.first_constraint >= AppParConstraint::TangencyPoint
            && self.last_constraint >= AppParConstraint::TangencyPoint
        {
            for j in 1..=ninc1 {
                the_a.set(length - 2 * self.ninc + j + 1, the_v1[j - 1]);
            }
            for j in 1..=self.ninc {
                the_a.set(length - self.ninc + j, the_v2[j - 1]);
            }
        } else if self.first_constraint >= AppParConstraint::TangencyPoint {
            for j in 1..=self.ninc {
                the_a.set(length - self.ninc + j, the_v1[j - 1]);
            }
        } else if self.last_constraint >= AppParConstraint::TangencyPoint {
            for j in 1..=self.ninc {
                the_a.set(length - self.ninc + j, the_v2[j - 1]);
            }
        }
    }

    /// OCCT Perform(Parameters, l1, l2) fixed-tangency resolution used by the
    /// ComputeCurve interpolation path (L794-1063): poles[2] = P1 + l1*V1t and
    /// poles[N-1] = P2 - l2*V2t are preset, the interior poles resinit..resfin
    /// (3..N-2) are solved by the reduced system.
    fn perform_tangency(&mut self) {
        let n_bcols = self.n_bcols;
        let nbpoles = self.nbpoles;
        let first_p = self.first_p;
        let last_p = self.last_p;
        // Preset the tangency poles (OCCT L820-836).
        if self.first_constraint >= AppParConstraint::TangencyPoint {
            for i in 0..n_bcols {
                self.mypoles[1][i] = self.mypoints[0][i] + self.lambda1 * self.vec1t[i];
            }
        }
        if self.last_constraint >= AppParConstraint::TangencyPoint {
            for i in 0..n_bcols {
                self.mypoles[nbpoles - 2][i] =
                    self.mypoints[self.mylastp - self.myfirstp][i] - self.lambda2 * self.vec2t[i];
            }
        }
        // B2 reduction: subtract the fixed poles (first, last, and the two
        // tangency poles).
        for j in first_p..=last_p {
            let dj = j - self.myfirstp;
            for i in 0..n_bcols {
                let mut v = self.mypoints[dj][i];
                if self.first_constraint >= AppParConstraint::PassPoint {
                    v -= self.a[dj][0] * self.mypoles[0][i];
                }
                if self.first_constraint >= AppParConstraint::TangencyPoint {
                    v -= self.a[dj][1] * self.mypoles[1][i];
                }
                if self.last_constraint >= AppParConstraint::PassPoint {
                    v -= self.a[dj][nbpoles - 1] * self.mypoles[nbpoles - 1][i];
                }
                if self.last_constraint >= AppParConstraint::TangencyPoint {
                    v -= self.a[dj][nbpoles - 2] * self.mypoles[nbpoles - 2][i];
                }
                self.b2[j - first_p][i] = v;
            }
        }
        let nincx = if self.resfin >= self.resinit {
            self.resfin - self.resinit + 1
        } else {
            0
        };
        if nincx < 1 {
            self.done = true;
            return;
        }
        let mut index = rcad_kernel::math::IntVec::new(nincx);
        self.search_index(&mut index);
        let mut mytab = rcad_kernel::math::MatD::new(nincx, n_bcols);
        let the_aa_len = index.get(nincx) as usize;
        let mut the_aa = rcad_kernel::math::VecD::new(the_aa_len);
        self.make_taa(&mut the_aa, &mut mytab);
        if rcad_kernel::math::lin::dactcl_decompose(&mut the_aa, &index.v, 1.0e-20) != 0 {
            self.done = false;
            return;
        }
        for j in 1..=n_bcols {
            let mut my_tabb = rcad_kernel::math::VecD::new(nincx);
            for (kk2, i) in (self.resinit..=self.resfin).enumerate() {
                my_tabb.set(kk2 + 1, mytab.get(i - self.resinit + 1, j));
            }
            if rcad_kernel::math::lin::dactcl_solve(&the_aa, &mut my_tabb, &index.v, 1.0e-20) != 0 {
                self.done = false;
                return;
            }
            let mut i2 = 0usize;
            for k in self.resinit..=self.resfin {
                self.mypoles[k - 1][j - 1] = my_tabb.get(i2 + 1);
                i2 += 1;
            }
        }
        self.done = true;
    }

    /// OCCT LeastSquare::SearchIndex (L1816-1859), Bezier case (no knots).
    fn search_index(&self, index: &mut rcad_kernel::math::IntVec) {
        let nincx = if self.resfin >= self.resinit {
            self.resfin - self.resinit + 1
        } else {
            0
        };
        index.set(1, 1);
        let mut l = 1usize;
        for i in 2..=nincx {
            l += 1;
            let v = index.get(l - 1) + i as i32;
            index.set(l, v);
        }
    }

    /// OCCT LeastSquare::MakeTAA(AA, myTAB) (L1705-1761), Bezier case: the
    /// reduced normal equations (triangle of the upper part) and the reduced
    /// right-hand sides.
    fn make_taa(&self, aa: &mut rcad_kernel::math::VecD, my_tab: &mut rcad_kernel::math::MatD) {
        let n_poles = if self.resfin >= self.resinit {
            self.resfin - self.resinit + 1
        } else {
            0
        };
        let mut the_a = rcad_kernel::math::MatD::new(n_poles, n_poles);
        for row in the_a.m.iter_mut() {
            for v in row.iter_mut() {
                *v = 0.0;
            }
        }
        for k in self.first_p..=self.last_p {
            let dk = k - self.myfirstp;
            let jinit = 1usize.max(self.resinit);
            let jfin = (self.deg + 1).min(self.resfin);
            for i in jinit..=jfin {
                let akj = self.a[dk][i - 1];
                for j in jinit..=i {
                    let v = the_a.get(i - self.resinit + 1, j - self.resinit + 1)
                        + self.a[dk][j - 1] * akj;
                    the_a.set(i - self.resinit + 1, j - self.resinit + 1, v);
                }
                for j in 1..=self.n_bcols {
                    let v = my_tab.get(i - self.resinit + 1, j) + akj * self.b2[k - self.first_p][j - 1];
                    my_tab.set(i - self.resinit + 1, j, v);
                }
            }
        }
        // Compress the triangle (Bezier: len = 2, single pass).
        let mut i2 = 1usize;
        let iinit = self.resinit;
        let mut jinit = self.resinit;
        let ifin = (self.deg + 1).min(self.resfin);
        for _k in 2..=2 {
            for i in iinit..=ifin {
                for j in jinit..=i {
                    let v = the_a.get(i - self.resinit + 1, j - self.resinit + 1);
                    aa.set(i2, v);
                    i2 += 1;
                }
            }
        }
    }

    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT BezierValue(): the resulting MultiCurve.
    pub fn bezier_value(&self) -> MultiCurve {
        let mut mc = MultiCurve::new(self.nbpoles, self.nb_p, self.nb_p2d);
        for (i, p) in self.mypoles.iter().enumerate() {
            let mut mp = MultiPoint::new(self.nb_p, self.nb_p2d);
            let mut i2 = 0;
            for j in 0..self.nb_p {
                mp.p3d[j] = DVec3::new(p[i2], p[i2 + 1], p[i2 + 2]);
                i2 += 3;
            }
            for j in 0..self.nb_p2d {
                mp.p2d[j] = DVec2::new(p[i2], p[i2 + 1]);
                i2 += 2;
            }
            mc.set_value(i + 1, mp);
        }
        mc
    }

    pub fn function_matrix(&self) -> &Vec<Vec<f64>> {
        &self.a
    }
    pub fn derivative_function_matrix(&self) -> &Vec<Vec<f64>> {
        &self.da
    }

    /// OCCT Error (L1214-1280): sum of squared distances and max errors.
    pub fn error(&self) -> (f64, f64, f64) {
        let mut f = 0.0;
        let mut max_e3d = 0.0f64;
        let mut max_e2d = 0.0f64;
        for r in self.first_p..=self.last_p {
            let data_idx = r - self.myfirstp;
            let mut k = 0usize;
            for _c in 0..self.nb_p {
                let mut fi = 0.0;
                for coord in 0..3 {
                    let mut aa = 0.0;
                    for p in 0..self.nbpoles {
                        aa += self.a[data_idx][p] * self.mypoles[p][k + coord];
                    }
                    let fx = aa - self.mypoints[data_idx][k + coord];
                    fi += fx * fx;
                }
                if fi > max_e3d {
                    max_e3d = fi;
                }
                f += fi;
                k += 3;
            }
            for _c in 0..self.nb_p2d {
                let mut fi = 0.0;
                for coord in 0..2 {
                    let mut aa = 0.0;
                    for p in 0..self.nbpoles {
                        aa += self.a[data_idx][p] * self.mypoles[p][k + coord];
                    }
                    let fx = aa - self.mypoints[data_idx][k + coord];
                    fi += fx * fx;
                }
                if fi > max_e2d {
                    max_e2d = fi;
                }
                f += fi;
                k += 2;
            }
        }
        (f, max_e3d.sqrt(), max_e2d.sqrt())
    }

    /// OCCT ErrorGradient (L1282-1369): sum of squared distances and the
    /// gradient with respect to the interior point parameters.
    pub fn error_gradient(&self, grad: &mut VecD, tol: f64) -> (f64, f64, f64) {
        let mut f = 0.0;
        let mut max_e3d = 0.0f64;
        let mut max_e2d = 0.0f64;
        for g in grad.v.iter_mut() {
            *g = 0.0;
        }
        let first_p = self.first_p;
        let last_p = self.last_p;
        for r in first_p..=last_p {
            let data_idx = r - self.myfirstp;
            let mut k = 0usize;
            let mut gr = 0.0;
            for _c in 0..self.nb_p {
                let mut fi = 0.0;
                for coord in 0..3 {
                    let mut aa = 0.0;
                    let mut daa = 0.0;
                    for p in 0..self.nbpoles {
                        aa += self.a[data_idx][p] * self.mypoles[p][k + coord];
                        daa += self.da[data_idx][p] * self.mypoles[p][k + coord];
                    }
                    let fx = aa - self.mypoints[data_idx][k + coord];
                    fi += fx * fx;
                    gr += 2.0 * daa * fx;
                }
                if fi > max_e3d {
                    max_e3d = fi;
                }
                f += fi;
                k += 3;
            }
            for _c in 0..self.nb_p2d {
                let mut fi = 0.0;
                for coord in 0..2 {
                    let mut aa = 0.0;
                    let mut daa = 0.0;
                    for p in 0..self.nbpoles {
                        aa += self.a[data_idx][p] * self.mypoles[p][k + coord];
                        daa += self.da[data_idx][p] * self.mypoles[p][k + coord];
                    }
                    let fx = aa - self.mypoints[data_idx][k + coord];
                    fi += fx * fx;
                    gr += 2.0 * daa * fx;
                }
                if fi > max_e2d {
                    max_e2d = fi;
                }
                f += fi;
                k += 2;
            }
            if gr.abs() <= tol {
                gr = 0.0;
            }
            grad.set(r - self.myfirstp + 1, gr);
        }
        (f, max_e3d.sqrt(), max_e2d.sqrt())
    }
}

/// Gaussian elimination with partial pivoting for a dense n x n system.
/// Returns None when the matrix is singular.  Semantic equivalent of the
/// OCCT math_Householder / DACTCL solves.
fn gauss_solve(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    let mut m = a.to_vec();
    let mut x = b.to_vec();
    for i in 0..n {
        // partial pivot
        let mut piv = i;
        for r in i + 1..n {
            if m[r][i].abs() > m[piv][i].abs() {
                piv = r;
            }
        }
        if m[piv][i].abs() < 1e-300 {
            return None;
        }
        m.swap(i, piv);
        x.swap(i, piv);
        for r in i + 1..n {
            let f = m[r][i] / m[i][i];
            if f == 0.0 {
                continue;
            }
            for c in i..n {
                m[r][c] -= f * m[i][c];
            }
            x[r] -= f * x[i];
        }
    }
    for i in (0..n).rev() {
        let mut s = x[i];
        for c in i + 1..n {
            s -= m[i][c] * x[c];
        }
        x[i] = s / m[i][i];
    }
    Some(x)
}

// ============================================================================
// AppParCurves_Function (AppParCurves_Function.gxx L36-696) — the function
// F = sum ||C(ui) - Ptli||^2 over the parameters, evaluated through the
// LeastSquare; used by the Gradient and the BFGS minimizer.
// ============================================================================

pub struct ParFunction<'a> {
    pub ml: &'a WLineAccess<'a>,
    pub my_parameters: VecD,
    pub my_multi_curve: MultiCurve,
    pub my_least_square: LeastSquare,
    pub f_val: f64,
    pub err3d: f64,
    pub err2d: f64,
    pub grad_val: VecD,
    pub first_p: usize,
    pub last_p: usize,
    pub adeb: usize,
    pub afin: usize,
    pub degre: usize,
    pub contraintes: bool,
    pub done: bool,
}

impl<'a> ParFunction<'a> {
    /// OCCT AppParCurves_Function constructor (L36-143).
    pub fn new(
        ml: &'a WLineAccess<'a>,
        first_point: usize,
        last_point: usize,
        constraints: &[ConstraintCouple],
        parameters: &VecD,
        deg: usize,
    ) -> Self {
        let mut adeb = first_point;
        let mut afin = last_point;
        let mut contraintes = false;
        for c in constraints {
            if c.index as usize == first_point {
                if c.constraint as i32 >= 1 {
                    adeb += 1;
                }
            } else if c.index as usize == last_point {
                if c.constraint as i32 >= 1 {
                    afin -= 1;
                }
            } else if c.constraint as i32 >= 1 {
                contraintes = true;
            }
        }
        let first_c = constraints
            .iter()
            .find(|c| c.index as usize == first_point)
            .map(|c| c.constraint)
            .unwrap_or(AppParConstraint::NoConstraint);
        let last_c = constraints
            .iter()
            .find(|c| c.index as usize == last_point)
            .map(|c| c.constraint)
            .unwrap_or(AppParConstraint::NoConstraint);
        let ls = LeastSquare::new(ml, first_point, last_point, first_c, last_c, deg + 1);
        let n = parameters.len();
        let mut my_parameters = VecD::new(n);
        for i in 1..=n {
            my_parameters.set(i, parameters.get(i));
        }
        if std::env::var("RCAD_TG_DEBUG").is_ok() {
            let pv: Vec<f64> = (1..=n.min(12)).map(|i| parameters.get(i)).collect();
            eprintln!("  [PAR] deg={} nbpoles={} first={} last={} pars={:?}", deg, deg + 1, first_point, last_point, pv);
        }
        ParFunction {
            ml,
            my_parameters,
            my_multi_curve: MultiCurve::new(deg + 1, ml.nb_p3d(), ml.nb_p2d()),
            my_least_square: ls,
            f_val: 0.0,
            err3d: 0.0,
            err2d: 0.0,
            grad_val: VecD::new(last_point - first_point + 1),
            first_p: first_point,
            last_p: last_point,
            adeb,
            afin,
            degre: deg,
            contraintes,
            done: false,
        }
    }

    /// OCCT Function::Value (L189-294).
    pub fn value(&mut self, x: &VecD) -> bool {
        for i in 1..=self.my_parameters.len() {
            self.my_parameters.set(i, x.get(i));
        }
        self.my_least_square.perform(&self.my_parameters);
        if !self.my_least_square.is_done() {
            self.done = false;
            return false;
        }
        if !self.contraintes {
            let (fv, e3, e2) = self.my_least_square.error();
            self.f_val = fv;
            self.err3d = e3;
            self.err2d = e2;
            self.done = true;
        } else {
            // Internal constraints (not reachable in the quadric WLine path);
            // evaluate F directly on the Bezier curve.
            self.my_multi_curve = self.my_least_square.bezier_value();
            let a = self.my_least_square.function_matrix().clone();
            let mut f = 0.0;
            let mut max_e3d = 0.0f64;
            let mut max_e2d = 0.0f64;
            let nb3d = self.ml.nb_p3d();
            let nb2d = self.ml.nb_p2d();
            for ci in 1..=nb3d {
                for i in self.adeb..=self.afin {
                    let mut aa = 0.0;
                    let mut bb = 0.0;
                    let mut cc = 0.0;
                    for j in 0..=self.degre {
                        let aij = a[i - self.first_p][j];
                        let pt = self.my_multi_curve.poles[j].point(ci);
                        aa += aij * pt.x;
                        bb += aij * pt.y;
                        cc += aij * pt.z;
                    }
                    let p = self.ml.value_p3d(i);
                    let fx = aa - p.x;
                    let fy = bb - p.y;
                    let fz = cc - p.z;
                    let fi = fx * fx + fy * fy + fz * fz;
                    if fi.sqrt() > max_e3d {
                        max_e3d = fi.sqrt();
                    }
                    f += fi;
                }
            }
            for ci in 1..=nb2d {
                for i in self.adeb..=self.afin {
                    let mut aa = 0.0;
                    let mut bb = 0.0;
                    for j in 0..=self.degre {
                        let aij = a[i - self.first_p][j];
                        let pt = self.my_multi_curve.poles[j].point2d(ci);
                        aa += aij * pt.x;
                        bb += aij * pt.y;
                    }
                    let p = self.ml.value_p2d(i)[ci - 1];
                    let fx = aa - p.x;
                    let fy = bb - p.y;
                    let fi = fx * fx + fy * fy;
                    if fi.sqrt() > max_e2d {
                        max_e2d = fi.sqrt();
                    }
                    f += fi;
                }
            }
            self.f_val = f;
            self.err3d = max_e3d;
            self.err2d = max_e2d;
            self.done = true;
        }
        true
    }

    /// OCCT Function::Perform (L296-646) — value + gradient.  The gradient
    /// with respect to the interior parameters (semantic equivalent).
    pub fn perform(&mut self, x: &VecD) {
        for i in 1..=self.my_parameters.len() {
            self.my_parameters.set(i, x.get(i));
        }
        self.my_least_square.perform(&self.my_parameters);
        if !self.my_least_square.is_done() {
            self.done = false;
            return;
        }
        for j in 1..=self.grad_val.len() {
            self.grad_val.set(j, 0.0);
        }
        if !self.contraintes {
            let (fv, e3, e2) = self.my_least_square.error_gradient(&mut self.grad_val, 0.0);
            self.f_val = fv;
            self.err3d = e3;
            self.err2d = e2;
        } else {
            // Constraint-free evaluation on the Bezier curve.
            let _ = self.value(x);
            self.f_val = self.f_val;
        }
        self.done = true;
    }

    pub fn curve_value(&mut self) -> MultiCurve {
        if !self.contraintes {
            self.my_multi_curve = self.my_least_square.bezier_value();
        }
        self.my_multi_curve.clone()
    }

    pub fn error(&self, i_point: usize, _curve: usize) -> f64 {
        // Recompute the per-point error from the least square matrix.
        let data_idx = i_point - self.my_least_square.myfirstp;
        let mut fi = 0.0;
        let mut k = 0usize;
        for _c in 0..self.my_least_square.nb_p {
            for coord in 0..3 {
                let mut aa = 0.0;
                for p in 0..self.my_least_square.nbpoles {
                    aa += self.my_least_square.a[data_idx][p]
                        * self.my_least_square.mypoles[p][k + coord];
                }
                let fx = aa - self.my_least_square.mypoints[data_idx][k + coord];
                fi += fx * fx;
            }
            k += 3;
        }
        for _c in 0..self.my_least_square.nb_p2d {
            for coord in 0..2 {
                let mut aa = 0.0;
                for p in 0..self.my_least_square.nbpoles {
                    aa += self.my_least_square.a[data_idx][p]
                        * self.my_least_square.mypoles[p][k + coord];
                }
                let fx = aa - self.my_least_square.mypoints[data_idx][k + coord];
                fi += fx * fx;
            }
            k += 2;
        }
        fi.sqrt()
    }

    pub fn max_error3d(&self) -> f64 {
        self.err3d
    }
    pub fn max_error2d(&self) -> f64 {
        self.err2d
    }
    pub fn new_parameters(&self) -> VecD {
        self.my_parameters.clone()
    }
}

// ============================================================================
// AppParCurves_Gradient (AppParCurves_Gradient.gxx L44-219) — the
// least-squares fit with the Rogers & Fog parameter projection and the
// optional BFGS refinement.
// ============================================================================

pub struct Gradient {
    pub scu: MultiCurve,
    pub done: bool,
    pub m_error3d: f64,
    pub m_error2d: f64,
    pub av_error: f64,
    pub par_error: VecD,
}

impl Gradient {
    /// OCCT AppParCurves_Gradient constructor (L44-219) with the BFGS
    /// refinement replaced by its semantic equivalent (a bounded number of
    /// Gauss-Newton steps on the interior parameters).
    pub fn new(
        ml: &WLineAccess,
        first_point: usize,
        last_point: usize,
        constraints: &[ConstraintCouple],
        parameters: &mut VecD,
        deg: usize,
        tol3d: f64,
        tol2d: f64,
        nb_iterations: i32,
    ) -> Self {
        let nb_p3d = ml.nb_p3d();
        let nb_p2d = ml.nb_p2d();
        let mut myf = ParFunction::new(ml, first_point, last_point, constraints, parameters, deg);
        let mut grad = Gradient {
            scu: MultiCurve::new(deg + 1, nb_p3d, nb_p2d),
            done: false,
            m_error3d: 0.0,
            m_error2d: 0.0,
            av_error: 0.0,
            par_error: VecD::new(last_point - first_point + 1),
        };
        let mut fval = 0.0;
        if !myf.value(parameters) {
            return grad;
        }
        fval = myf.f_val;
        grad.scu = myf.curve_value();
        // OCCT AppParCurves_Gradient.gxx L99-125: storage of curve poles for
        // projection, converted to power-basis coefficients via
        // BSplCLib::PolesCoefficients (TheCoef / TheCoef2d).
        let mut the_coef: Vec<DVec3> = Vec::with_capacity((deg + 1) * nb_p3d);
        let mut tab_pole = Vec::with_capacity(deg + 1);
        for k in 0..nb_p3d {
            grad.scu.curve(k + 1, &mut tab_pole);
            let tab_coef = bezier_poles_to_coeffs(&tab_pole);
            the_coef.extend_from_slice(&tab_coef);
        }
        let mut the_coef2d: Vec<DVec2> = Vec::with_capacity((deg + 1) * nb_p2d);
        let mut tab_pole2d = Vec::with_capacity(deg + 1);
        for k in 0..nb_p2d {
            grad.scu.curve2d(nb_p3d + k + 1, &mut tab_pole2d);
            let tab_coef2d = bezier_poles_to_coeffs_2d(&tab_pole2d);
            the_coef2d.extend_from_slice(&tab_coef2d);
        }
        // OCCT L131-175: Rogers & Fog projection iteration (no D2 needed).
        for j in first_point + 1..last_point {
            let uf = parameters.get(j - first_point + 1);
            let p3 = ml.value_p3d(j);
            let p2 = ml.value_p2d(j);
            let mut fu = 0.0;
            let mut dfu = 0.0;
            let mut i2 = 0usize;
            for k in 0..nb_p3d {
                // OCCT L148-151: TabCoef(l) = TheCoef(l + i2);
                // BSplCLib::CoefsD1(UF, TabCoef, NoWeights, Pt, V1) =
                // CacheD1 (BSplCLib_CurveComputation.pxx L1307-1346) with
                // SpanLenght = 1: a Horner power-basis evaluation.
                let tab_coef = &the_coef[i2..i2 + deg + 1];
                i2 += deg + 1;
                let mut pt = DVec3::ZERO;
                let mut v1 = DVec3::ZERO;
                for d in 0..3 {
                    let cx: Vec<f64> = tab_coef.iter().map(|c| c[d]).collect();
                    let (v, dv) = eval_polynomial_d1(&cx, uf);
                    pt[d] = v;
                    v1[d] = dv;
                }
                // OCCT L152-154: MyV = gp_Vec(Pt, TabP(k)); FU += MyV * V1;
                // DFU += V1.SquareMagnitude();
                let myv = p3 - pt;
                fu += myv.dot(v1);
                dfu += v1.length_squared();
            }
            let mut i2 = 0usize;
            for k in 0..nb_p2d {
                let tab_coef = &the_coef2d[i2..i2 + deg + 1];
                i2 += deg + 1;
                let mut pt2 = DVec2::ZERO;
                let mut v12 = DVec2::ZERO;
                for d in 0..2 {
                    let cx: Vec<f64> = tab_coef.iter().map(|c| c[d]).collect();
                    let (v, dv) = eval_polynomial_d1(&cx, uf);
                    pt2[d] = v;
                    v12[d] = dv;
                }
                // OCCT L163-165: MyV2d = gp_Vec2d(Pt2d, TabP2d(k));
                // FU += MyV2d * V12d; DFU += V12d.SquareMagnitude();
                let myv2 = p2[k] - pt2;
                fu += myv2.dot(v12);
                dfu += v12.length_squared();
            }
            // OCCT L168-174: DFU >= RealEpsilon(); DU clamped to +/-5e-02.
            if dfu >= 2.220446049250313e-16 {
                let mut du = fu / dfu;
                du = du.abs().min(5.0e-02).copysign(du);
                let uf2 = uf + du;
                parameters.set(j - first_point + 1, uf2);
            }
        }
        if !myf.value(parameters) {
            grad.scu = MultiCurve::new(deg + 1, nb_p3d, nb_p2d);
            grad.done = false;
            return grad;
        }
        fval = myf.f_val;
        grad.m_error3d = myf.max_error3d();
        grad.m_error2d = myf.max_error2d();
        if std::env::var("RCAD_TG_DEBUG").is_ok() {
            let pv: Vec<f64> = (1..=parameters.len().min(12)).map(|i| parameters.get(i)).collect();
            eprintln!("  [GRAD] deg={} first={} last={} err3d={:.3e} err2d={:.3e} after_proj pars={:?}", deg, first_point, last_point, grad.m_error3d, grad.m_error2d, pv);
        }
        if grad.m_error3d <= tol3d && grad.m_error2d <= tol2d {
            grad.done = true;
            grad.scu = myf.curve_value();
        } else if nb_iterations != 0 {
            // OCCT L191-199: AppParCurves_Gradient_BFGS — math_BFGS on the
            // interior parameters (math_BFGS.cxx L327-443) with the analytic
            // value + gradient of the ParFunction.
            let n = parameters.len();
            let start: Vec<f64> = (1..=n).map(|i| parameters.get(i)).collect();
            let mut myf_ref = &mut myf;
            let opt = rcad_kernel::math::opt::bfgs_minimize_occt(
                n,
                &start,
                1.0e-7,
                nb_iterations,
                1.0e-7,
                &mut |x: &[f64]| -> Option<(f64, Vec<f64>)> {
                    let mut xd = rcad_kernel::math::VecD::new(n);
                    for i in 0..n {
                        xd.set(i + 1, x[i]);
                    }
                    myf_ref.perform(&xd);
                    if myf_ref.done {
                        let g: Vec<f64> = (1..=myf_ref.grad_val.len())
                            .map(|i| myf_ref.grad_val.get(i))
                            .collect();
                        Some((myf_ref.f_val, g))
                    } else {
                        None
                    }
                },
            );
            if let Some(best) = opt {
                for i in 0..n {
                    parameters.set(i + 1, best[i]);
                }
            }
            grad.m_error3d = myf.max_error3d();
            grad.m_error2d = myf.max_error2d();
            grad.scu = myf.curve_value();
        }
        // Average error (OCCT L201-211).
        let mut av = 0.0;
        for j in first_point..=last_point {
            let mut pe: f64 = 0.0;
            let n_curves = nb_p3d + nb_p2d;
            for k in 1..=n_curves {
                pe = pe.max(myf.error(j, k));
            }
            grad.par_error.set(j - first_point + 1, pe);
            av += pe;
        }
        grad.av_error = av / (last_point - first_point + 1) as f64;
        grad.m_error3d = myf.max_error3d();
        grad.m_error2d = myf.max_error2d();
        if grad.m_error3d <= tol3d && grad.m_error2d <= tol2d {
            grad.done = true;
        }
        grad
    }

    pub fn value(&self) -> MultiCurve {
        self.scu.clone()
    }
    pub fn is_done(&self) -> bool {
        self.done
    }
    pub fn max_error3d(&self) -> f64 {
        self.m_error3d
    }
    pub fn max_error2d(&self) -> f64 {
        self.m_error2d
    }
}

// ============================================================================
// ApproxInt_KnotTools — 1:1 translation of ApproxInt_KnotTools.cxx
// ============================================================================

/// OCCT ApproxInt_KnotTools::EvalCurv (L36-84): curvature of an n-dim curve
/// from the first and second derivative vectors V1, V2 (length `dim`).
fn eval_curv(dim: usize, v1: &[f64], v2: &[f64]) -> f64 {
    let mut mp = 0.0;
    for i in 1..dim {
        for j in 0..i {
            let p = v1[i] * v2[j] - v1[j] * v2[i];
            mp += p * p;
        }
    }
    let mut q = 0.0;
    for i in 0..dim {
        q += v1[i] * v1[i];
    }
    if q < 1.0 / 2.0e100 {
        // Singularity: imitate a curvature jump (OCCT L62-75).
        return 2.0e100;
    }
    let q = q.min(2.0e100);
    let q = q * q * q;
    (mp / q).sqrt()
}

/// OCCT ApproxInt_KnotTools::BuildCurvature (L88-167): discrete curvature of
/// the n-dim curve `theCoords` (per-point dim coordinates) at each parameter.
fn build_curvature(coords: &[f64], dim: usize, n: usize, pars: &[f64]) -> (Vec<f64>, f64) {
    let mut curv = vec![0.0; n];
    let mut max_curv = 0.0;
    if n < 3 {
        return (curv, max_curv);
    }
    // First point: Lagrange through points 0,1,2 evaluated at par[0].
    {
        let mut val = vec![0.0; 3 * dim];
        for m in 0..dim {
            val[m] = coords[m];
            val[dim + m] = coords[dim + m];
            val[2 * dim + m] = coords[2 * dim + m];
        }
        let (_, d1, d2) = eval_lagrange(&val, dim, pars[0], pars[1], pars[2], pars[0]);
        curv[0] = eval_curv(dim, &d1, &d2);
        max_curv = max_curv.max(curv[0]);
    }
    // Interior points: Lagrange through i-1, i, i+1 evaluated at par[i].
    for i in 1..n - 1 {
        let mut val = vec![0.0; 3 * dim];
        for m in 0..dim {
            val[m] = coords[(i - 1) * dim + m];
            val[dim + m] = coords[i * dim + m];
            val[2 * dim + m] = coords[(i + 1) * dim + m];
        }
        let (_, d1, d2) = eval_lagrange(&val, dim, pars[i - 1], pars[i], pars[i + 1], pars[i]);
        curv[i] = eval_curv(dim, &d1, &d2);
        max_curv = max_curv.max(curv[i]);
    }
    // Last point: Lagrange through n-3, n-2, n-1 evaluated at par[n-1].
    {
        let mut val = vec![0.0; 3 * dim];
        for m in 0..dim {
            val[m] = coords[(n - 3) * dim + m];
            val[dim + m] = coords[(n - 2) * dim + m];
            val[2 * dim + m] = coords[(n - 1) * dim + m];
        }
        let (_, d1, d2) = eval_lagrange(
            &val,
            dim,
            pars[n - 3],
            pars[n - 2],
            pars[n - 1],
            pars[n - 1],
        );
        curv[n - 1] = eval_curv(dim, &d1, &d2);
        max_curv = max_curv.max(curv[n - 1]);
    }
    (curv, max_curv)
}

/// OCCT PLib::EvalLagrange(Parameter, DerivativeRequest=2, Degree=2, Dim,
/// Values, Parameters, Results) (PLib.cxx L1122-1249): Newton divided-difference
/// evaluation of the degree-2 Lagrange interpolant through the three nodes
/// (t0, t1, t2) with the values `val`, evaluated at the parameter t.  Returns
/// (value, first derivative, second derivative), each of length `dim`.
fn eval_lagrange(val: &[f64], dim: usize, t0: f64, t1: f64, t2: f64, t: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let degree = 2usize;
    let local_request = 2usize; // min(DerivativeRequest, Degree)
    let par = [t0, t1, t2];
    // Build the divided differences array.
    let mut dd = vec![0.0; (degree + 1) * dim];
    for i in 0..(degree + 1) * dim {
        dd[i] = val[i];
    }
    let mut ok = true;
    'outer: for ii in (0..=degree).rev() {
        for jj in ((degree - ii + 1)..=degree).rev() {
            let index = jj * dim;
            let index1 = index - dim;
            for kk in 0..dim {
                dd[index + kk] -= dd[index1 + kk];
            }
            let difference = par[jj] - par[jj + ii - degree - 1];
            if difference.abs() < 2.2250738585072014e-308 {
                // OCCT: |difference| < RealSmall() -> ReturnCode = 1; goto FINISH.
                ok = false;
                break 'outer;
            }
            let difference = 1.0 / difference;
            for kk in 0..dim {
                dd[index + kk] *= difference;
            }
        }
    }
    // Evaluate the divided-difference polynomial (Newton form).
    let mut res = vec![0.0; (local_request + 1) * dim];
    if ok {
        let index = degree * dim;
        for kk in 0..dim {
            res[kk] = dd[index + kk];
        }
        for i in dim..(local_request + 1) * dim {
            res[i] = 0.0;
        }
        for ii in (1..=degree).rev() {
            let difference = t - par[ii - 1];
            for jj in (1..=local_request).rev() {
                let index = jj * dim;
                let index1 = index - dim;
                for kk in 0..dim {
                    res[index + kk] *= difference;
                    res[index + kk] += res[index1 + kk] * jj as f64;
                }
            }
            let index = (ii - 1) * dim;
            for kk in 0..dim {
                res[kk] *= difference;
                res[kk] += dd[index + kk];
            }
        }
    }
    let res0 = res[0..dim].to_vec();
    let res1 = res[dim..2 * dim].to_vec();
    let res2 = res[2 * dim..3 * dim].to_vec();
    (res0, res1, res2)
}

/// OCCT ApproxInt_KnotTools::InsKnotBefI (L452-552): try to insert a knot
/// between theInds(theI-1) and theInds(theI) (1-based sequence access);
/// `chk_curv` selects the curvature-change or the angular criteria.
/// Returns the inserted index (mid) or None.
fn ins_knot_bef_i(
    the_i: usize,
    curv: &[f64],
    coords: &[f64],
    dim: usize,
    inds: &[usize],
    chk_curv: bool,
) -> Option<usize> {
    let an_ind1 = inds[the_i]; // theInds(theI)
    let an_ind = inds[the_i - 1]; // theInds(theI-1)
    if an_ind1 - an_ind == 1 {
        return None;
    }
    let a_limit_curvature_change = 3.0;
    let a_sin_coeff2 = 0.09549150281252627; // (3 - sqrt(5)) / 8
    // OCCT L466-467: curv = 0.5 * (theCurv(anInd) + theCurv(anInd1)).
    let curv_half = 0.5 * (curv[an_ind] + curv[an_ind1]);
    let mut mid = 0usize;
    let mut j = an_ind + 1;
    while j < an_ind1 {
        mid = 0;
        // I: curvature change criteria.
        if curv[j] > 1e-7 && curv[an_ind] > 1e-7 {
            if curv[j] / curv[an_ind] > a_limit_curvature_change
                || curv[j] / curv[an_ind] < 1.0 / a_limit_curvature_change
            {
                mid = j;
                inds_insert(&mut inds.to_vec(), the_i, mid);
                return Some(mid);
            }
        }
        // II: angular criteria.
        let ac = curv[j - 1];
        let ac1 = curv[j];
        if (curv_half >= ac && curv_half <= ac1) || (curv_half >= ac1 && curv_half <= ac) {
            if (curv_half - ac).abs() < (curv_half - ac1).abs() {
                mid = j - 1;
            } else {
                mid = j;
            }
        }
        if mid == an_ind {
            mid += 1;
        }
        if mid == an_ind1 {
            mid -= 1;
        }
        if mid > 0 {
            if chk_curv {
                let ici = an_ind * dim;
                let ici1 = an_ind1 * dim;
                let icm = mid * dim;
                let mut v1 = vec![0.0; dim];
                let mut v2 = vec![0.0; dim];
                let mut m1 = 0.0;
                let mut m2 = 0.0;
                let mut mp = 0.0;
                for i in 0..dim {
                    v1[i] = coords[icm + i] - coords[ici + i];
                    m1 += v1[i] * v1[i];
                    v2[i] = coords[ici1 + i] - coords[icm + i];
                    m2 += v2[i] * v2[i];
                }
                for i in 1..dim {
                    for jj in 0..i {
                        let p = v1[i] * v2[jj] - v1[jj] * v2[i];
                        mp += p * p;
                    }
                }
                if mp > a_sin_coeff2 * m1 * m2 {
                    inds_insert(&mut inds.to_vec(), the_i, mid);
                    return Some(mid);
                }
            } else {
                inds_insert(&mut inds.to_vec(), the_i, mid);
                return Some(mid);
            }
        }
        j += 1;
    }
    None
}

fn inds_insert(inds: &mut Vec<usize>, pos: usize, val: usize) {
    inds.insert(pos, val);
}

/// OCCT ApproxInt_KnotTools::ComputeKnotInds (L171-327).
fn compute_knot_inds(coords: &[f64], dim: usize, n: usize, pars: &[f64]) -> (Vec<usize>, Vec<usize>) {
    // I: create the discrete curvature.
    let (a_curv, a_max_curv) = build_curvature(coords, dim, n, pars);
    let mut inds: Vec<usize> = Vec::new();
    let mut feature_inds: Vec<usize> = Vec::new();
    inds.push(0);
    if a_max_curv <= 1e-7 {
        // Linear case.
        inds.push(n - 1);
        return (inds, feature_inds);
    }
    // II: find extremas of curvature.
    let eps = 1.0e-9;
    let eps1 = 1.0e3 * eps;
    for i in 1..n - 1 {
        let d1 = a_curv[i] - a_curv[i - 1];
        let d2 = a_curv[i] - a_curv[i + 1];
        let ad1 = d1.abs();
        let ad2 = d2.abs();
        if d1 * d2 > 0.0 && ad1 > eps && ad2 > eps {
            if *inds.last().unwrap() != i {
                inds.push(i);
                feature_inds.push(i);
            }
        } else if (ad1 < eps && ad2 > eps1) || (ad1 > eps1 && ad2 < eps) {
            if *inds.last().unwrap() != i {
                inds.push(i);
                feature_inds.push(i);
            }
        }
    }
    if n - 1 != *inds.last().unwrap() {
        inds.push(n - 1);
    }
    // III: put knots in monotone intervals of curvature (OCCT L241-252).
    // OCCT: i = 1; do { i++; Ok = InsKnotBefI(i, ...); if (Ok) i--; } while (i < theInds.Length());
    // Here `i` is the 1-based sequence position theI; ins_knot_bef_i takes the 0-based the_i = theI - 1.
    let mut i = 1usize;
    loop {
        i += 1;
        if let Some(mid) = ins_knot_bef_i(i - 1, &a_curv, coords, dim, &inds, true) {
            inds.insert(i - 1, mid);
            i -= 1;
        }
        if i >= inds.len() {
            break; // while (i < theInds.Length())
        }
    }
    // IV: checking feature points (OCCT L254-325).
    // OCCT: j = 2; for (; j <= theInds.Length() - 1;) — j is a 1-based sequence
    // position.  rcad `j` is the 0-based position, so it starts at 1 and runs
    // while j < inds.len() - 1.
    let mut j = 1usize;
    let mut fi = 0usize;
    while fi < feature_inds.len() {
        let an_ind = feature_inds[fi];
        let mut inserted = false;
        while j < inds.len() - 1 {
            if inds[j] == an_ind {
                let an_ind_prev = inds[j - 1];
                let an_ind_next = inds[j + 1];
                let ici = an_ind_prev * dim;
                let ici1 = an_ind_next * dim;
                let icm = an_ind * dim;
                let mut v1 = vec![0.0; dim];
                let mut v2 = vec![0.0; dim];
                let mut m1 = 0.0;
                let mut m2 = 0.0;
                let mut mp = 0.0;
                for k in 0..dim {
                    v1[k] = coords[icm + k] - coords[ici + k];
                    m1 += v1[k] * v1[k];
                    v2[k] = coords[ici1 + k] - coords[icm + k];
                    m2 += v2[k] * v2[k];
                }
                for k in 1..dim {
                    for l in 0..k {
                        let p = v1[k] * v2[l] - v1[l] * v2[k];
                        mp += p * p;
                    }
                }
                if mp > 0.09549150281252627 * m1 * m2 {
                    let d1 = (a_curv[an_ind] - a_curv[an_ind_prev]).abs();
                    let d2 = (a_curv[an_ind] - a_curv[an_ind_next]).abs();
                    if d1 > d2 {
                        if let Some(mid) = ins_knot_bef_i(j, &a_curv, coords, dim, &inds, false) {
                            inds.insert(j, mid);
                            inserted = true;
                            j += 1;
                        } else {
                            break;
                        }
                    } else if let Some(mid) = ins_knot_bef_i(j + 1, &a_curv, coords, dim, &inds, false)
                    {
                        inds.insert(j + 1, mid);
                        inserted = true;
                    } else {
                        break;
                    }
                } else {
                    j += 1;
                    break;
                }
                if inserted {
                    // re-check this feature with the updated sequence
                    continue;
                }
            } else {
                j += 1;
            }
        }
        fi += 1;
    }
    (inds, feature_inds)
}

/// OCCT ApproxInt_KnotTools::FilterKnots (L331-448).
fn filter_knots(inds: &mut Vec<usize>, min_nb_pnts: usize) -> Vec<usize> {
    let a_max_nb_pnts = 15 * min_nb_pnts;
    let a_min_nb_step = min_nb_pnts / 2;
    // I: filter too big number of points per knot interval.
    let mut i = 1usize;
    while i < inds.len() {
        let nbint = inds[i] - inds[i - 1] + 1;
        if nbint <= a_max_nb_pnts {
            i += 1;
            continue;
        } else {
            let ind = inds[i - 1] + nbint / 2;
            inds.insert(i, ind);
        }
    }
    // II: filter points with too small amount of points per knot interval.
    let mut lknots: Vec<usize> = Vec::new();
    i = 1;
    lknots.push(inds[0]);
    let mut an_inds_prev = inds[0];
    i = 2;
    while i <= inds.len() {
        if inds[i - 1] - an_inds_prev <= min_nb_pnts {
            if i != inds.len() {
                let mut an_idx = i + 1;
                while an_idx <= inds.len() {
                    if inds[an_idx - 1] - an_inds_prev >= min_nb_pnts {
                        break;
                    }
                    an_idx += 1;
                }
                an_idx -= 1;
                let a_mid_idx = (inds[an_idx - 1] + an_inds_prev) / 2;
                if (a_mid_idx as i64 - an_inds_prev as i64) < min_nb_pnts as i64
                    && (a_mid_idx as i64 - inds[an_idx - 1] as i64) < min_nb_pnts as i64
                    && inds[an_idx - 1] - an_inds_prev >= a_min_nb_step
                {
                    if inds[an_idx - 1] - an_inds_prev > 2 * min_nb_pnts {
                        lknots.push(an_inds_prev + min_nb_pnts);
                        an_inds_prev = an_inds_prev + min_nb_pnts;
                        i = an_idx - 1;
                    } else {
                        if inds[an_idx - 2] - an_inds_prev >= min_nb_pnts / 2 {
                            lknots.push(inds[an_idx - 2]);
                            an_inds_prev = inds[an_idx - 2];
                            i = an_idx - 1;
                            if inds[an_idx - 1] - inds[an_idx - 2] <= min_nb_pnts / 2 {
                                *lknots.last_mut().unwrap() = inds[an_idx - 1];
                                an_inds_prev = inds[an_idx - 1];
                                i = an_idx;
                            }
                        } else {
                            lknots.push(inds[an_idx - 1]);
                            an_inds_prev = inds[an_idx - 1];
                            i = an_idx;
                        }
                    }
                } else if an_idx == inds.len() && lknots.len() >= 2 {
                    let a_last_good_idx = lknots[lknots.len() - 2];
                    if inds[inds.len() - 1] - 2 * min_nb_pnts >= a_last_good_idx {
                        *lknots.last_mut().unwrap() = inds[inds.len() - 1] - min_nb_pnts;
                        lknots.push(inds[inds.len() - 1]);
                        an_inds_prev = inds[an_idx - 1];
                        i = an_idx;
                    }
                }
            }
            i += 1;
            continue;
        } else {
            lknots.push(inds[i - 1]);
            an_inds_prev = inds[i - 1];
            i += 1;
        }
    }
    // III: fill the last knot.
    if lknots.len() < 2 {
        lknots.push(*inds.last().unwrap());
    } else {
        if *lknots.last().unwrap() < *inds.last().unwrap() {
            *lknots.last_mut().unwrap() = *inds.last().unwrap();
        }
    }
    lknots
}

/// OCCT ApproxInt_KnotTools::BuildKnots (L556-640) — the full knot sequence
/// builder on the part points.  `coords` holds per-point concatenated
/// coordinates (dim per point), `pars` the normalized parameters (unused by
/// the knot criteria themselves; kept for signature parity), `min_nb_pnts`
/// the OCCT aMinNbPnts (= myNbPntMax).  Returns the knot POINT INDICES.
pub fn build_knots(
    coords: &[f64],
    dim: usize,
    pars: &[f64],
    min_nb_pnts: usize,
) -> Vec<usize> {
    let n = pars.len();
    if dim == 0 || n < 2 {
        return Vec::new();
    }
    let (mut draft, _feat) = compute_knot_inds(coords, dim, n, pars);
    filter_knots(&mut draft, min_nb_pnts)
}

/// OCCT ApproxInt_KnotTools::MaxParamRatio (L644-665).
fn max_param_ratio(pars: &[f64]) -> f64 {
    let mut a_max_ratio: f64 = 0.0;
    for i in 1..pars.len() - 1 {
        let a_denom = pars[i] - pars[i - 1];
        if a_denom.abs() < 2.220446049250313e-16 {
            // OCCT L652: std::abs(aDenom) < Precision::Computational() (RealEpsilon()).
            continue;
        }
        let mut a_rat = (pars[i + 1] - pars[i]) / a_denom;
        if a_rat > 0.0 && a_rat < 1.0 {
            a_rat = 1.0 / a_rat;
        }
        a_max_ratio = a_max_ratio.max(a_rat);
    }
    a_max_ratio
}

/// OCCT ApproxInt_Approx::Parameters (ApproxInt_Approx.gxx L88-162) —
/// chord-length / centripetal / iso-parametric parameterization of the
/// MultiLine points on [firstP, lastP], normalized to [0, 1].
pub fn approx_parameters(
    ml: &WLineAccess,
    first_p: usize,
    last_p: usize,
    par: ApproxParamType,
) -> VecD {
    let n = last_p - first_p + 1;
    let mut the_parameters = VecD::new(n);
    if par == ApproxParamType::ChordLength || par == ApproxParamType::Centripetal {
        let nb_p3d = ml.nb_p3d();
        let nb_p2d = ml.nb_p2d();
        let mynb_p3d = if nb_p3d == 0 { 1 } else { nb_p3d };
        let mynb_p2d = if nb_p2d == 0 { 1 } else { nb_p2d };
        the_parameters.set(1, 0.0);
        for i in first_p + 1..=last_p {
            let p1 = ml.value_p3d(i - 1);
            let p2 = ml.value_p3d(i);
            let p2d1 = ml.value_p2d(i - 1);
            let p2d2 = ml.value_p2d(i);
            let mut dist = 0.0;
            for j in 0..nb_p3d {
                let _ = (mynb_p3d, j);
                let a = p1;
                let b = p2;
                dist += (b - a).length_squared();
            }
            for j in 0..nb_p2d {
                let _ = (mynb_p2d, j);
                let a = p2d1[j];
                let b = p2d2[j];
                dist += (b - a).length_squared();
            }
            dist = dist.sqrt();
            if par == ApproxParamType::ChordLength {
                the_parameters.set(i - first_p + 1, the_parameters.get(i - first_p) + dist);
            } else {
                the_parameters.set(i - first_p + 1, the_parameters.get(i - first_p) + dist.sqrt());
            }
        }
        let last = the_parameters.get(n);
        if last > 0.0 {
            for i in 1..=n {
                the_parameters.set(i, the_parameters.get(i) / last);
            }
        }
    } else {
        for i in 0..n {
            the_parameters.set(i + 1, i as f64 / (n - 1) as f64);
        }
    }
    the_parameters
}

/// OCCT ApproxInt_KnotTools::DefineParType (L669-848): choose the
/// parameterization by the curvature profile of the part.
pub fn define_par_type(
    ml: &WLineAccess,
    the_fpar: usize,
    the_lpar: usize,
    the_approx_xyz: bool,
    the_approx_u1v1: bool,
    the_approx_u2v2: bool,
) -> ApproxParamType {
    if the_lpar - the_fpar == 1 {
        return ApproxParamType::IsoParametric;
    }
    // Build the concatenated per-point coordinates.
    let a_dim = (if the_approx_xyz { 3 } else { 0 })
        + (if the_approx_u1v1 { 2 } else { 0 })
        + (if the_approx_u2v2 { 2 } else { 0 });
    let a_length = the_lpar - the_fpar + 1;
    let mut a_coords = vec![0.0f64; a_length * a_dim];
    for i in 0..a_length {
        let idx = the_fpar + i;
        let mut j = i * a_dim;
        let p3 = ml.value_p3d(idx);
        let p2 = ml.value_p2d(idx);
        if the_approx_xyz {
            a_coords[j] = p3.x;
            j += 1;
            a_coords[j] = p3.y;
            j += 1;
            a_coords[j] = p3.z;
            j += 1;
        }
        if the_approx_u1v1 {
            a_coords[j] = p2[0].x;
            j += 1;
            a_coords[j] = p2[0].y;
            j += 1;
        }
        if the_approx_u2v2 {
            let k = if p2.len() >= 2 { 1 } else { 0 };
            a_coords[j] = p2[k].x;
            j += 1;
            a_coords[j] = p2[k].y;
            j += 1;
        }
    }
    // Analysis of curvature.
    let a_crit_rat = 500.0;
    let a_crit_par_rat = 100.0;
    let a_pars = approx_parameters(ml, the_fpar, the_lpar, ApproxParamType::ChordLength);
    let mut a_par_type = ApproxParamType::ChordLength;
    let (a_curv, a_max_curv) = build_curvature(&a_coords, a_dim, a_length, &a_pars.v);
    if a_max_curv < 1.0e-9 || a_max_curv >= 1.0e100 {
        // OCCT L800: aMaxCurv < Precision::PConfusion() || Precision::IsPositiveInfinite(aMaxCurv).
        return a_par_type;
    }
    let mut a_mid_curv = 0.0;
    let eps = 1.0f64 * 2.220446049250313e-16;
    let mut j = 0usize;
    for i in 0..a_length {
        if a_max_curv - a_curv[i] < eps {
            continue;
        }
        j += 1;
        a_mid_curv += a_curv[i];
    }
    if j > 1 {
        a_mid_curv /= j as f64;
    }
    if a_mid_curv <= eps {
        return a_par_type;
    }
    let a_rat = a_max_curv / a_mid_curv;
    if a_rat > a_crit_rat {
        if a_rat > 5.0 * a_crit_rat {
            a_par_type = ApproxParamType::Centripetal;
        } else {
            let a_par_rat = max_param_ratio(&a_pars.v);
            if a_par_rat > a_crit_par_rat {
                a_par_type = ApproxParamType::Centripetal;
            }
        }
    }
    a_par_type
}

