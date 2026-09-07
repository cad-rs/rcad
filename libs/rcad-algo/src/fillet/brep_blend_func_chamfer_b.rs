//! OCCT BlendFunc_ChAsym (TKFillet/BlendFunc) — 1:1 port of
//! BlendFunc_ChAsym.hxx (L34-222) + BlendFunc_ChAsym.cxx (whole file
//! L33-768); and BlendFunc_ChAsymInv — 1:1 port of BlendFunc_ChAsymInv.hxx
//! (L26-88) + BlendFunc_ChAsymInv.cxx (whole file L25-403).  The structs are
//! defined in [`super::brep_blend_func_chamfer`]; this file carries their
//! method bodies and trait wiring.
//!
//! Architecture mappings: `class BlendFunc_ChAsym : public Blend_Function`
//! is expressed by implementing the [`BlendFunction`] trait and the
//! `math_FunctionSetWithDerivatives` base; `class BlendFunc_ChAsymInv :
//! public Blend_FuncInv` by the [`BlendFuncInv`] trait.  `tcurv` (OCCT
//! handle aliasing `curv` until Set(First, Last)) maps to an owned trimmed
//! [`Curve3`] copy plus an accessor.

use glam::{DVec2, DVec3};

use rcad_kernel::core::precision::is_infinite_value;
use rcad_kernel::geom::{Curve3, CurveEval as _, Curve2d, Curve2dEval as _, Surface3, SurfaceEval as _, TrimmedCurve3};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;
use rcad_kernel::math::gp::Lin;
use rcad_kernel::math::math_gauss::MathGauss;
use rcad_kernel::math::{GeomAbsShape, MatD, VecD};

use super::brep_blend_func::blend_func_next_shape;
use super::brep_blend_func_chamfer::{BlendFuncChAsym, BlendFuncChAsymInv};
use super::brep_blend_func_inv::BlendFuncInv;
use super::brep_blend_function::{BlendAppFunction, BlendFunction};
use super::brep_blend_point::BlendPoint;

impl<'a> BlendFuncChAsym<'a> {
    /// Architecture mapping: OCCT `tcurv` is a handle aliasing `curv` until
    /// Set(First, Last) replaces it with a trimmed copy (ChAsym.cxx L74).
    #[inline]
    fn tcurv(&self) -> &Curve3 {
        self.tcurv.as_ref().unwrap_or(self.curv)
    }

    /// OCCT BlendFunc_ChAsym(S1, S2, C) (BlendFunc_ChAsym.cxx L33-50).
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, c: &'a Curve3) -> Self {
        BlendFuncChAsym {
            surf1: s1,
            surf2: s2,
            curv: c,
            tcurv: None, // OCCT: tcurv(C) — alias, see tcurv().
            param: 0.0,
            dist1: f64::MAX, // OCCT: RealLast()
            angle: f64::MAX,
            tgang: f64::MAX,
            nplan: DVec3::ZERO,
            pt1: DVec3::ZERO,
            tsurf1: DVec3::ZERO,
            pt2: DVec3::ZERO,
            fx: [0.0; 4],
            dx: vec![vec![0.0; 4]; 4],
            istangent: true,
            tg1: DVec3::ZERO,
            tg12d: DVec2::ZERO,
            tg2: DVec3::ZERO,
            tg22d: DVec2::ZERO,
            choix: 0,
            distmin: f64::MAX,
        }
    }

    /// OCCT NbEquations() (BlendFunc_ChAsym.cxx L54-57) — returns 4.
    pub fn nb_equations(&self) -> usize {
        4
    }

    /// OCCT Set(Param) (BlendFunc_ChAsym.cxx L61-64).
    pub fn set_param(&mut self, param: f64) {
        self.param = param;
    }

    /// OCCT Set(First, Last) (BlendFunc_ChAsym.cxx L72-75) — segmente la
    /// courbe a sa partie utile.  La precision est prise arbitrairement
    /// petite !?
    pub fn set_interval(&mut self, first: f64, last: f64) {
        // OCCT: tcurv = curv->Trim(First, Last, 1.e-12);
        self.tcurv = Some(Curve3::Trimmed(TrimmedCurve3 {
            curve: Box::new(self.curv.clone()),
            first,
            last,
        }));
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BlendFunc_ChAsym.cxx L79-85).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        tolerance[0] = self.surf1.u_resolution(tol);
        tolerance[1] = self.surf1.v_resolution(tol);
        tolerance[2] = self.surf2.u_resolution(tol);
        tolerance[3] = self.surf2.v_resolution(tol);
    }

    /// OCCT GetBounds(InfBound, SupBound) (BlendFunc_ChAsym.cxx L89-109).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        inf_bound[0] = self.surf1.default_domain()[0]; // FirstUParameter
        inf_bound[1] = self.surf1.default_domain()[2]; // FirstVParameter
        inf_bound[2] = self.surf2.default_domain()[0];
        inf_bound[3] = self.surf2.default_domain()[2];
        sup_bound[0] = self.surf1.default_domain()[1]; // LastUParameter
        sup_bound[1] = self.surf1.default_domain()[3]; // LastVParameter
        sup_bound[2] = self.surf2.default_domain()[1];
        sup_bound[3] = self.surf2.default_domain()[3];

        for i in 0..4 {
            if !is_infinite_value(inf_bound[i]) && !is_infinite_value(sup_bound[i]) {
                let range = sup_bound[i] - inf_bound[i];
                inf_bound[i] -= range;
                sup_bound[i] += range;
            }
        }
    }

    /// OCCT IsSolution(Sol, Tol) (BlendFunc_ChAsym.cxx L113-206).
    #[allow(unreachable_code)] // the math_SVD fallback is a pending kernel gap
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT: math_Vector valsol(1, 4), secmember(1, 4);
        //       math_Matrix gradsol(1, 4, 1, 4);
        let mut valsol = [0.0f64; 4];
        let mut secmember = [0.0f64; 4];
        let mut gradsol = vec![vec![0.0f64; 4]; 4];

        // OCCT: tcurv->D2(param, ptgui, d1gui, d2gui);
        let ptgui = self.tcurv().point_at(self.param);
        let d1gui = self.tcurv().derivative_at(self.param);
        let d2gui = self.tcurv().derivative2_at(self.param);
        let mut normg = d1gui.length();
        let mut np = d1gui.normalize_or_zero();
        // OCCT: dnp = (d2gui - np.Dot(d2gui) * np) / Normg;
        let mut dnp = (d2gui - np.dot(d2gui) * np) / normg;

        if self.choix % 2 != 0 {
            np = -np;
            dnp = -dnp;
            normg = -normg;
        }

        // OCCT: surf1->D1(Sol(1), Sol(2), pt1, d1u1, d1v1);
        let (p1, d1u1, d1v1) = {
            let (p, du, dv) = self.surf1.derivatives(sol[0], sol[1]);
            (p, du, dv)
        };
        self.pt1 = p1;
        let nsurf1 = d1u1.cross(d1v1);
        self.tsurf1 = nsurf1.cross(np);
        // OCCT: dwtsurf1 = Nsurf1.Crossed(dnp);
        let dwtsurf1 = nsurf1.cross(dnp);

        // OCCT: surf2->D1(Sol(3), Sol(4), pt2, d1u2, d1v2);
        let (p2, d1u2, d1v2) = {
            let (p, du, dv) = self.surf2.derivatives(sol[2], sol[3]);
            (p, du, dv)
        };
        self.pt2 = p2;

        // OCCT: gp_Vec pguis1(ptgui, pt1), pguis2(ptgui, pt2);
        let pguis1 = self.pt1 - ptgui;
        let pguis2 = self.pt2 - ptgui;
        let s1s2 = self.pt2 - self.pt1;
        let psca_inv = 1.0 / self.tsurf1.dot(s1s2);
        let nordu1 = d1u1.length();
        let nordv1 = d1v1.length();

        // OCCT: temp = 2. * (Nordu1 + Nordv1) * s1s2.Magnitude() + 2. * Nordu1 * Nordv1;
        let temp = 2.0 * (nordu1 + nordv1) * s1s2.length() + 2.0 * nordu1 * nordv1;

        self.values(sol, &mut valsol, &mut gradsol);

        if valsol[0].abs() < tol
            && valsol[1].abs() < tol
            && valsol[2].abs() < 2.0 * self.dist1 * tol
            && valsol[3].abs() < tol * (1.0 + self.tgang) * psca_inv.abs() * temp
        {
            secmember[0] = normg - dnp.dot(pguis1);
            secmember[1] = normg - dnp.dot(pguis2);
            secmember[2] = -2.0 * d1gui.dot(pguis1);

            // OCCT: CrossVec = tsurf1.Crossed(s1s2);
            let cross_vec = self.tsurf1.cross(s1s2);
            let f4 = np.dot(cross_vec) * psca_inv;
            let mut temp = dnp.dot(cross_vec) + np.dot(dwtsurf1.cross(s1s2));

            temp -= f4 * dwtsurf1.dot(s1s2);
            secmember[3] = psca_inv * temp;

            // OCCT: math_Gauss Resol(gradsol, maxpiv);
            let mut a = MatD::new(4, 4);
            for r in 1..=4 {
                for c in 1..=4 {
                    a.set(r, c, gradsol[r - 1][c - 1]);
                }
            }
            let resol = MathGauss::new(&a);

            if resol.is_done() {
                let mut x = VecD::new(4);
                for i in 1..=4 {
                    x.set(i, secmember[i - 1]);
                }
                resol.solve(&mut x);
                for i in 1..=4 {
                    secmember[i - 1] = x.get(i);
                }
                self.istangent = false;
            } else {
                // OCCT L177-188: math_SVD SingRS(gradsol); if
                // (SingRS.IsDone()) { SingRS.Solve(DEDT, secmember, 1.e-6); }.
                // GAP (plan 0.6): rcad-kernel exposes no public math_SVD yet
                // (only private svd helpers in function_set_root.rs); the SVD
                // fallback is pending kernel support and reports a tangency
                // point like the OCCT !IsDone() path.
                self.istangent = true;
            }

            if !self.istangent {
                // OCCT: tg1.SetLinearForm(secmember(1), d1u1, secmember(2), d1v1);
                self.tg1 = secmember[0] * d1u1 + secmember[1] * d1v1;
                self.tg2 = secmember[2] * d1u2 + secmember[3] * d1v2;
                self.tg12d = DVec2::new(secmember[0], secmember[1]);
                self.tg22d = DVec2::new(secmember[2], secmember[3]);
            }

            self.distmin = self.distmin.min(self.pt1.distance(self.pt2));

            return true;
        }

        self.istangent = true;
        false
    }

    /// OCCT GetMinimalDistance() (BlendFunc_ChAsym.cxx L210-213).
    pub fn get_minimal_distance(&self) -> f64 {
        self.distmin
    }

    /// OCCT ComputeValues(X, DegF, DegL) (BlendFunc_ChAsym.cxx L217-306).
    pub fn compute_values(&mut self, x: &[f64], deg_f: i32, deg_l: i32) -> bool {
        if deg_f > deg_l {
            return false;
        }

        // OCCT: tcurv->D1(param, ptgui, d1gui);
        let ptgui = self.tcurv().point_at(self.param);
        let d1gui = self.tcurv().derivative_at(self.param);
        self.nplan = d1gui.normalize_or_zero();
        let mut np = self.nplan;

        if self.choix % 2 != 0 {
            np = -np;
        }

        let d1u1: DVec3;
        let d1v1: DVec3;
        let d2u1: DVec3;
        let d2v1: DVec3;
        let d2uv1: DVec3;
        let d1u2: DVec3;
        let d1v2: DVec3;

        if (deg_f == 0) && (deg_l == 0) {
            // OCCT: surf1->D1(X(1), X(2), pt1, d1u1, d1v1);
            let (p1, du1, dv1) = self.surf1.derivatives(x[0], x[1]);
            self.pt1 = p1;
            d1u1 = du1;
            d1v1 = dv1;
            d2u1 = DVec3::ZERO;
            d2v1 = DVec3::ZERO;
            d2uv1 = DVec3::ZERO;
            // OCCT: pt2 = surf2->Value(X(3), X(4));
            self.pt2 = self.surf2.point_at(x[2], x[3]);
            d1u2 = DVec3::ZERO;
            d1v2 = DVec3::ZERO;
        } else {
            // OCCT: surf1->D2(X(1), X(2), pt1, d1u1, d1v1, d2u1, d2v1, d2uv1);
            let (p1, du1, dv1, d2u1v, d2uv1v, d2v1v) = self.surf1.derivatives2(x[0], x[1]);
            self.pt1 = p1;
            d1u1 = du1;
            d1v1 = dv1;
            d2u1 = d2u1v;
            d2uv1 = d2uv1v;
            d2v1 = d2v1v;
            // OCCT: surf2->D1(X(3), X(4), pt2, d1u2, d1v2);
            let (p2, du2, dv2) = self.surf2.derivatives(x[2], x[3]);
            self.pt2 = p2;
            d1u2 = du2;
            d1v2 = dv2;
        }

        let nsurf1 = d1u1.cross(d1v1);
        self.tsurf1 = nsurf1.cross(np);

        // OCCT: gp_Vec nps1(ptgui, pt1), s1s2(pt1, pt2);
        let nps1 = self.pt1 - ptgui;
        let s1s2 = self.pt1 - self.pt2;
        let psca_inv = 1.0 / self.tsurf1.dot(s1s2);
        let f4 = np.dot(self.tsurf1.cross(s1s2)) * psca_inv;

        if deg_f == 0 {
            // OCCT: Dist = ptgui.XYZ().Dot(np.XYZ());
            let dist = ptgui.dot(np);

            self.fx[0] = self.pt1.dot(np) - dist;
            self.fx[1] = self.pt2.dot(np) - dist;
            self.fx[2] = self.dist1 * self.dist1 - nps1.length_squared();
            self.fx[3] = self.tgang - f4;
        }

        if deg_l == 1 {
            // OCCT: d1utsurf1 = (d2u1.Crossed(d1v1) + d1u1.Crossed(d2uv1)).Crossed(np);
            let d1utsurf1 = (d2u1.cross(d1v1) + d1u1.cross(d2uv1)).cross(np);
            let d1vtsurf1 = (d2uv1.cross(d1v1) + d1u1.cross(d2v1)).cross(np);

            self.dx[0][0] = np.dot(d1u1);
            self.dx[0][1] = np.dot(d1v1);
            self.dx[0][2] = 0.0;
            self.dx[0][3] = 0.0;

            self.dx[1][0] = 0.0;
            self.dx[1][1] = 0.0;
            self.dx[1][2] = np.dot(d1u2);
            self.dx[1][3] = np.dot(d1v2);

            // OCCT: tempVec = -2. * nps1;
            let temp_vec = -2.0 * nps1;
            self.dx[2][0] = d1u1.dot(temp_vec);
            self.dx[2][1] = d1v1.dot(temp_vec);
            self.dx[2][2] = 0.0;
            self.dx[2][3] = 0.0;

            let mut temp = f4 * (d1utsurf1.dot(s1s2) - self.tsurf1.dot(d1u1));
            temp += np.dot(self.tsurf1.cross(d1u1) - d1utsurf1.cross(s1s2));
            self.dx[3][0] = temp * psca_inv;

            let mut temp = f4 * (d1vtsurf1.dot(s1s2) - self.tsurf1.dot(d1v1));
            temp += np.dot(self.tsurf1.cross(d1v1) - d1vtsurf1.cross(s1s2));
            self.dx[3][1] = temp * psca_inv;

            let temp = f4 * self.tsurf1.dot(d1u2) - np.dot(self.tsurf1.cross(d1u2));
            self.dx[3][2] = temp * psca_inv;

            let temp = f4 * self.tsurf1.dot(d1v2) - np.dot(self.tsurf1.cross(d1v2));
            self.dx[3][3] = temp * psca_inv;
        }

        true
    }

    /// OCCT Value(X, F) (BlendFunc_ChAsym.cxx L310-315).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let error = self.compute_values(x, 0, 0);
        f.copy_from_slice(&self.fx);
        error
    }

    /// OCCT Derivatives(X, D) (BlendFunc_ChAsym.cxx L319-324).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        let error = self.compute_values(x, 1, 1);
        // OCCT: D = DX;
        for (drow, xrow) in d.iter_mut().zip(self.dx.iter()) {
            drow.clone_from(xrow);
        }
        error
    }

    /// OCCT Values(X, F, D) (BlendFunc_ChAsym.cxx L328-334).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        let error = self.compute_values(x, 0, 1);
        f.copy_from_slice(&self.fx);
        // OCCT: D = DX;
        for (drow, xrow) in d.iter_mut().zip(self.dx.iter()) {
            drow.clone_from(xrow);
        }
        error
    }

    /// OCCT PointOnS1() (BlendFunc_ChAsym.cxx L338-341).
    pub fn point_on_s1(&self) -> DVec3 {
        self.pt1
    }

    /// OCCT PointOnS2() (BlendFunc_ChAsym.cxx L345-348).
    pub fn point_on_s2(&self) -> DVec3 {
        self.pt2
    }

    /// OCCT IsTangencyPoint() (BlendFunc_ChAsym.cxx L352-355).
    pub fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// OCCT TangentOnS1() (BlendFunc_ChAsym.cxx L359-366).
    pub fn tangent_on_s1(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ChAsym::TangentOnS1");
        }
        self.tg1
    }

    /// OCCT Tangent2dOnS1() (BlendFunc_ChAsym.cxx L370-377).
    pub fn tangent_2d_on_s1(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ChAsym::Tangent2dOnS1");
        }
        self.tg12d
    }

    /// OCCT TangentOnS2() (BlendFunc_ChAsym.cxx L381-388).
    pub fn tangent_on_s2(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ChAsym::TangentOnS2");
        }
        self.tg2
    }

    /// OCCT Tangent2dOnS2() (BlendFunc_ChAsym.cxx L392-399).
    pub fn tangent_2d_on_s2(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ChAsym::Tangent2dOnS2");
        }
        self.tg22d
    }

    /// OCCT TwistOnS1() (BlendFunc_ChAsym.cxx L403-410).
    pub fn twist_on_s1(&self) -> bool {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ChAsym::TwistOnS1");
        }
        self.tg1.dot(self.nplan) < 0.0
    }

    /// OCCT TwistOnS2() (BlendFunc_ChAsym.cxx L414-421).
    pub fn twist_on_s2(&self) -> bool {
        if self.istangent {
            panic!("Standard_DomainError: BlendFunc_ChAsym::TwistOnS2");
        }
        self.tg2.dot(self.nplan) < 0.0
    }

    /// OCCT Tangent(U1, V1, U2, V2, TgFirst, TgLast, NormFirst, NormLast)
    /// (BlendFunc_ChAsym.cxx L429-480) — TgF, NmF et TgL, NmL les tangentes
    /// et normales respectives aux surfaces S1 et S2.
    #[allow(clippy::too_many_arguments)]
    pub fn tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_f: &mut DVec3,
        tg_l: &mut DVec3,
        nm_f: &mut DVec3,
        nm_l: &mut DVec3,
    ) {
        let mut rev_f = false;
        let mut rev_l = false;

        // OCCT: tcurv->D1(param, ptgui, d1gui); np = d1gui.Normalized();
        let ptgui = self.tcurv().point_at(self.param);
        let d1gui = self.tcurv().derivative_at(self.param);
        let np = d1gui.normalize_or_zero();
        let _ = ptgui;

        // OCCT: surf1->D1(U1, V1, Pt1, d1u1, d1v1); NmF = d1u1.Crossed(d1v1);
        let (_, d1u1, d1v1) = self.surf1.derivatives(u1, v1);
        *nm_f = d1u1.cross(d1v1);

        // OCCT: surf2->D1(U2, V2, Pt2, d1u2, d1v2); NmL = d1u2.Crossed(d1v2);
        let (_, d1u2, d1v2) = self.surf2.derivatives(u2, v2);
        *nm_l = d1u2.cross(d1v2);

        *tg_f = np.cross(*nm_f).normalize_or_zero();
        *tg_l = np.cross(*nm_l).normalize_or_zero();

        if (self.choix == 2) || (self.choix == 5) {
            rev_f = true;
            rev_l = true;
        }

        if (self.choix == 4) || (self.choix == 7) {
            rev_l = true;
        }

        if (self.choix == 3) || (self.choix == 8) {
            rev_f = true;
        }

        if rev_f {
            *tg_f = -*tg_f;
        }
        if rev_l {
            *tg_l = -*tg_l;
        }
    }

    /// OCCT Section(Param, U1, V1, U2, V2, Pdeb, Pfin, C)
    /// (BlendFunc_ChAsym.cxx L484-502) — utile pour une visu rapide et
    /// approximative de la surface.
    #[allow(clippy::too_many_arguments)]
    pub fn section(
        &mut self,
        _param: f64,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        pdeb: &mut f64,
        pfin: &mut f64,
        c: &mut Lin,
    ) {
        let pt1 = self.surf1.point_at(u1, v1);
        let pt2 = self.surf2.point_at(u2, v2);
        // OCCT: const gp_Dir dir(gp_Vec(Pt1, Pt2));
        let dir = (pt2 - pt1).normalize_or_zero();

        c.pos = pt1;
        c.dir = dir;

        *pdeb = 0.0;
        // OCCT: Pfin = ElCLib::Parameter(C, Pt2);
        *pfin = (pt2 - c.pos).dot(c.dir);
    }

    /// OCCT IsRational() (BlendFunc_ChAsym.cxx L506-509).
    pub fn is_rational(&self) -> bool {
        false
    }

    /// OCCT GetSectionSize() (BlendFunc_ChAsym.cxx L513-516).
    pub fn get_section_size(&self) -> f64 {
        panic!("Standard_NotImplemented: BlendFunc_ChAsym::GetSectionSize()");
    }

    /// OCCT GetMinimalWeight(Weights) (BlendFunc_ChAsym.cxx L520-523).
    pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
        for w in weigths.iter_mut() {
            *w = 1.0;
        }
    }

    /// OCCT NbIntervals(S) (BlendFunc_ChAsym.cxx L527-530).
    pub fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        self.curv.nb_intervals(blend_func_next_shape(s))
    }

    /// OCCT Intervals(T, S) (BlendFunc_ChAsym.cxx L534-537).
    pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let mut intervals = Vec::new();
        self.curv.intervals(&mut intervals, blend_func_next_shape(s));
        for (dst, src) in t.iter_mut().zip(intervals) {
            *dst = src;
        }
    }

    /// OCCT GetShape(NbPoles, NbKnots, Degree, NbPoles2d)
    /// (BlendFunc_ChAsym.cxx L541-547).
    pub fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        *nb_poles = 2;
        *nb_poles_2d = 2;
        *nb_knots = 2;
        *degree = 1;
    }

    /// OCCT GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1D)
    /// (BlendFunc_ChAsym.cxx L553-560).
    pub fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        _surf_tol: f64,
        _angle_tol: f64,
        tol3d: &mut [f64],
        _tol1d: &mut [f64],
    ) {
        for v in tol3d.iter_mut() {
            *v = bound_tol;
        }
    }

    /// OCCT Knots(TKnots) (BlendFunc_ChAsym.cxx L564-568).
    pub fn knots(&mut self, tknots: &mut [f64]) {
        tknots[0] = 0.0;
        tknots[1] = 1.0;
    }

    /// OCCT Mults(TMults) (BlendFunc_ChAsym.cxx L572-576).
    pub fn mults(&mut self, tmults: &mut [i32]) {
        tmults[0] = 2;
        tmults[1] = 2;
    }

    /// OCCT Section(P, Poles, Poles2d, Weights)
    /// (BlendFunc_ChAsym.cxx L580-605).
    pub fn section_simple(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        weigths: &mut [f64],
    ) {
        let prm = p.parameter();
        let low = 0usize; // OCCT: Poles.Lower()
        let upp = poles.len() - 1; // OCCT: Poles.Upper()
        let mut x = [0.0f64; 4];
        let mut f = [0.0f64; 4];

        let (u1, v1) = p.parameters_on_s1();
        let (u2, v2) = p.parameters_on_s2();
        x[0] = u1;
        x[1] = v1;
        x[2] = u2;
        x[3] = v2;
        poles_2d[0] = DVec2::new(u1, v1);
        poles_2d[poles_2d.len() - 1] = DVec2::new(u2, v2);

        self.set_param(prm);
        self.value(&x, &mut f);
        poles[low] = self.point_on_s1();
        poles[upp] = self.point_on_s2();
        weigths[low] = 1.0;
        weigths[upp] = 1.0;
    }

    /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weights, DWeights)
    /// (BlendFunc_ChAsym.cxx L609-723) — used for the first and last section.
    #[allow(clippy::too_many_arguments)]
    #[allow(unreachable_code)] // the math_SVD fallback is a pending kernel gap
    pub fn section_d1(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
    ) -> bool {
        // OCCT: math_Vector Sol(1, 4), valsol(1, 4), secmember(1, 4);
        //       math_Matrix gradsol(1, 4, 1, 4);
        let mut sol = [0.0f64; 4];
        let mut valsol = [0.0f64; 4];
        let mut secmember = [0.0f64; 4];
        let mut gradsol = vec![vec![0.0f64; 4]; 4];
        let prm = p.parameter();
        let low = 0usize; // OCCT: Poles.Lower()
        let upp = poles.len() - 1; // OCCT: Poles.Upper()

        let (u1, v1) = p.parameters_on_s1();
        sol[0] = u1;
        sol[1] = v1;
        let (u2, v2) = p.parameters_on_s2();
        sol[2] = u2;
        sol[3] = v2;
        self.set_param(prm);

        poles_2d[0] = DVec2::new(sol[0], sol[1]);
        poles_2d[poles_2d.len() - 1] = DVec2::new(sol[2], sol[3]);
        poles[low] = self.point_on_s1();
        poles[upp] = self.point_on_s2();
        weigths[low] = 1.0;
        weigths[upp] = 1.0;

        // OCCT: tcurv->D2(param, ptgui, d1gui, d2gui);
        let ptgui = self.tcurv().point_at(self.param);
        let d1gui = self.tcurv().derivative_at(self.param);
        let d2gui = self.tcurv().derivative2_at(self.param);
        let mut normg = d1gui.length();
        let mut np = d1gui.normalize_or_zero();
        let mut dnp = (d2gui - np.dot(d2gui) * np) / normg;

        if self.choix % 2 != 0 {
            np = -np;
            dnp = -dnp;
            normg = -normg;
        }

        // OCCT: surf1->D1(Sol(1), Sol(2), pt1, d1u1, d1v1);
        let (p1, d1u1, d1v1) = {
            let (pp, du, dv) = self.surf1.derivatives(sol[0], sol[1]);
            (pp, du, dv)
        };
        self.pt1 = p1;
        let nsurf1 = d1u1.cross(d1v1);
        self.tsurf1 = nsurf1.cross(np);
        let dwtsurf1 = nsurf1.cross(dnp);

        // OCCT: surf2->D1(Sol(3), Sol(4), pt2, d1u2, d1v2);
        let (p2, d1u2, d1v2) = {
            let (pp, du, dv) = self.surf2.derivatives(sol[2], sol[3]);
            (pp, du, dv)
        };
        self.pt2 = p2;

        let pguis1 = self.pt1 - ptgui;
        let pguis2 = self.pt2 - ptgui;
        let s1s2 = self.pt2 - self.pt1;
        let psca_inv = 1.0 / self.tsurf1.dot(s1s2);
        let nordu1 = d1u1.length();
        let nordv1 = d1v1.length();

        // OCCT L664: temp is assigned here and overwritten at L674 before any
        // use — translated literally.
        let _temp = 2.0 * (nordu1 + nordv1) * s1s2.length() + 2.0 * nordu1 * nordv1;

        self.values(&sol, &mut valsol, &mut gradsol);

        secmember[0] = normg - dnp.dot(pguis1);
        secmember[1] = normg - dnp.dot(pguis2);
        secmember[2] = -2.0 * d1gui.dot(pguis1);

        let cross_vec = self.tsurf1.cross(s1s2);
        let f4 = np.dot(cross_vec) * psca_inv;
        let mut temp = dnp.dot(cross_vec) + np.dot(dwtsurf1.cross(s1s2));
        temp -= f4 * dwtsurf1.dot(s1s2);
        secmember[3] = psca_inv * temp;

        // OCCT: math_Gauss Resol(gradsol, maxpiv);
        let mut a = MatD::new(4, 4);
        for r in 1..=4 {
            for c in 1..=4 {
                a.set(r, c, gradsol[r - 1][c - 1]);
            }
        }
        let resol = MathGauss::new(&a);

        if resol.is_done() {
            let mut x = VecD::new(4);
            for i in 1..=4 {
                x.set(i, secmember[i - 1]);
            }
            resol.solve(&mut x);
            for i in 1..=4 {
                secmember[i - 1] = x.get(i);
            }
            self.istangent = false;
        } else {
            // OCCT L687-698: math_SVD fallback — see IsSolution for the GAP
            // note (plan 0.6, math_SVD pending in rcad-kernel).
            self.istangent = true;
        }

        if !self.istangent {
            self.tg1 = secmember[0] * d1u1 + secmember[1] * d1v1;
            self.tg2 = secmember[2] * d1u2 + secmember[3] * d1v2;
            self.tg12d = DVec2::new(secmember[0], secmember[1]);
            self.tg22d = DVec2::new(secmember[2], secmember[3]);
        }

        self.distmin = self.distmin.min(self.pt1.distance(self.pt2));

        if !self.istangent {
            let t2d1 = self.tangent_2d_on_s1();
            d_poles_2d[0] = DVec2::new(t2d1.x, t2d1.y);
            let t2d2 = self.tangent_2d_on_s2();
            d_poles_2d[poles_2d.len() - 1] = DVec2::new(t2d2.x, t2d2.y);

            d_poles[low] = self.tangent_on_s1();
            d_poles[upp] = self.tangent_on_s2();
            d_weigths[low] = 0.0;
            d_weigths[upp] = 0.0;
        }

        !self.istangent
    }

    /// OCCT Section(P, Poles, DPoles, D2Poles, Poles2d, DPoles2d, D2Poles2d,
    /// Weights, DWeights, D2Weights) (BlendFunc_ChAsym.cxx L727-739).
    #[allow(clippy::too_many_arguments)]
    pub fn section_d2(
        &mut self,
        _p: &BlendPoint,
        _poles: &mut [DVec3],
        _d_poles: &mut [DVec3],
        _d2_poles: &mut [DVec3],
        _poles_2d: &mut [DVec2],
        _d_poles_2d: &mut [DVec2],
        _d2_poles_2d: &mut [DVec2],
        _weigths: &mut [f64],
        _d_weigths: &mut [f64],
        _d2_weigths: &mut [f64],
    ) -> bool {
        false
    }

    /// OCCT Resolution(IC2d, Tol, TolU, TolV) (BlendFunc_ChAsym.cxx
    /// L743-758).
    pub fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        if ic_2d == 1 {
            *tol_u = self.surf1.u_resolution(tol);
            *tol_v = self.surf1.v_resolution(tol);
        } else {
            *tol_u = self.surf2.u_resolution(tol);
            *tol_v = self.surf2.v_resolution(tol);
        }
    }

    /// OCCT Set(Dist1, Angle, Choix) (BlendFunc_ChAsym.cxx L762-768) — sets
    /// the distances and the angle.
    pub fn set(&mut self, dist1: f64, angle: f64, choix: i32) {
        self.dist1 = dist1.abs();
        self.angle = angle;
        self.tgang = angle.tan();
        self.choix = choix;
    }
}

impl<'a> FunctionSetWithDerivatives for BlendFuncChAsym<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: inherited from Blend_Function::NbVariables (returns 4).
        BlendFunction::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendFuncChAsym::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendFuncChAsym::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncChAsym::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncChAsym::values(self, x, f, df)
    }
}

impl<'a> BlendAppFunction for BlendFuncChAsym<'a> {
    fn set_param(&mut self, param: f64) {
        BlendFuncChAsym::set_param(self, param)
    }

    fn set_interval(&mut self, first: f64, last: f64) {
        BlendFuncChAsym::set_interval(self, first, last)
    }

    // OCCT Blend_Function.cxx L24-32 — Pnt1/Pnt2 delegate to
    // PointOnS1/PointOnS2.
    fn pnt1(&self) -> DVec3 {
        BlendFunction::pnt1(self)
    }

    fn pnt2(&self) -> DVec3 {
        BlendFunction::pnt2(self)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendFuncChAsym::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendFuncChAsym::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendFuncChAsym::is_solution(self, sol, tol)
    }

    fn get_minimal_distance(&self) -> f64 {
        BlendFuncChAsym::get_minimal_distance(self)
    }

    fn is_rational(&self) -> bool {
        BlendFuncChAsym::is_rational(self)
    }

    fn get_section_size(&self) -> f64 {
        BlendFuncChAsym::get_section_size(self)
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        BlendFuncChAsym::get_minimal_weight(self, weigths)
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        BlendFuncChAsym::nb_intervals(self, s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        BlendFuncChAsym::intervals(self, t, s)
    }

    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        BlendFuncChAsym::get_shape(self, nb_poles, nb_knots, degree, nb_poles_2d)
    }

    fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        BlendFuncChAsym::get_approx_tolerance(self, bound_tol, surf_tol, angle_tol, tol3d, tol1d)
    }

    fn knots(&mut self, tknots: &mut [f64]) {
        BlendFuncChAsym::knots(self, tknots)
    }

    fn mults(&mut self, tmults: &mut [i32]) {
        BlendFuncChAsym::mults(self, tmults)
    }

    fn section_d1(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
    ) -> bool {
        BlendFuncChAsym::section_d1(
            self, p, poles, d_poles, poles_2d, d_poles_2d, weigths, d_weigths,
        )
    }

    fn section(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        weigths: &mut [f64],
    ) {
        BlendFuncChAsym::section_simple(self, p, poles, poles_2d, weigths)
    }

    fn section_d2(
        &mut self,
        p: &BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        d2_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        d2_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
        d2_weigths: &mut [f64],
    ) -> bool {
        BlendFuncChAsym::section_d2(
            self,
            p,
            poles,
            d_poles,
            d2_poles,
            poles_2d,
            d_poles_2d,
            d2_poles_2d,
            weigths,
            d_weigths,
            d2_weigths,
        )
    }

    fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        BlendFuncChAsym::resolution(self, ic_2d, tol, tol_u, tol_v)
    }
}

impl<'a> BlendFunction for BlendFuncChAsym<'a> {
    fn point_on_s1(&self) -> DVec3 {
        BlendFuncChAsym::point_on_s1(self)
    }

    fn point_on_s2(&self) -> DVec3 {
        BlendFuncChAsym::point_on_s2(self)
    }

    fn is_tangency_point(&self) -> bool {
        BlendFuncChAsym::is_tangency_point(self)
    }

    fn tangent_on_s1(&self) -> DVec3 {
        BlendFuncChAsym::tangent_on_s1(self)
    }

    fn tangent_2d_on_s1(&self) -> DVec2 {
        BlendFuncChAsym::tangent_2d_on_s1(self)
    }

    fn tangent_on_s2(&self) -> DVec3 {
        BlendFuncChAsym::tangent_on_s2(self)
    }

    fn tangent_2d_on_s2(&self) -> DVec2 {
        BlendFuncChAsym::tangent_2d_on_s2(self)
    }

    fn twist_on_s1(&self) -> bool {
        BlendFuncChAsym::twist_on_s1(self)
    }

    fn twist_on_s2(&self) -> bool {
        BlendFuncChAsym::twist_on_s2(self)
    }

    fn tangent(
        &self,
        u1: f64,
        v1: f64,
        u2: f64,
        v2: f64,
        tg_first: &mut DVec3,
        tg_last: &mut DVec3,
        norm_first: &mut DVec3,
        norm_last: &mut DVec3,
    ) {
        BlendFuncChAsym::tangent(self, u1, v1, u2, v2, tg_first, tg_last, norm_first, norm_last)
    }
}

impl<'a> BlendFuncChAsymInv<'a> {
    /// OCCT BlendFunc_ChAsymInv(S1, S2, C) (BlendFunc_ChAsymInv.cxx L25-39).
    pub fn new(s1: &'a Surface3, s2: &'a Surface3, c: &'a Curve3) -> Self {
        BlendFuncChAsymInv {
            surf1: s1,
            surf2: s2,
            dist1: f64::MAX, // OCCT: RealLast()
            angle: f64::MAX,
            tgang: f64::MAX,
            curv: c,
            csurf: None,
            choix: 0,
            first: false,
            fx: [0.0; 4],
            dx: vec![vec![0.0; 4]; 4],
        }
    }

    /// OCCT Set(Dist1, Angle, Choix) (BlendFunc_ChAsymInv.cxx L43-49).
    pub fn set(&mut self, dist1: f64, angle: f64, choix: i32) {
        self.dist1 = dist1.abs();
        self.angle = angle;
        self.tgang = angle.tan();
        self.choix = choix;
    }

    /// OCCT NbEquations() (BlendFunc_ChAsymInv.cxx L53-56) — returns 4.
    pub fn nb_equations(&self) -> usize {
        4
    }

    /// OCCT Set(OnFirst, C) (BlendFunc_ChAsymInv.cxx L60-64).
    pub fn set_curve_on_surface(&mut self, on_first: bool, c: &'a Curve2d) {
        self.first = on_first;
        self.csurf = Some(c);
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BlendFunc_ChAsymInv.cxx L68-82).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        // OCCT L70: Tolerance(1) = csurf->Resolution(Tol).
        // GAP (plan 0.6): rcad-kernel has no Adaptor2d_Curve2d::Resolution
        // equivalent yet; the call panics with the pending marker until the
        // kernel exposes it.
        tolerance[0] = super::brep_blend_func_chamfer::adaptor2d_curve2d_resolution_pending();
        // OCCT L71: Tolerance(2) = curv->Resolution(Tol).
        tolerance[1] = self.curv.resolution(tol);
        if self.first {
            tolerance[2] = self.surf2.u_resolution(tol);
            tolerance[3] = self.surf2.v_resolution(tol);
        } else {
            tolerance[2] = self.surf1.u_resolution(tol);
            tolerance[3] = self.surf1.v_resolution(tol);
        }
    }

    /// OCCT GetBounds(InfBound, SupBound) (BlendFunc_ChAsymInv.cxx L86-131).
    pub fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        let csurf = self.csurf.expect("csurf");
        inf_bound[0] = csurf.default_domain()[0]; // FirstParameter
        inf_bound[1] = self.curv.default_domain()[0];
        sup_bound[0] = csurf.default_domain()[1]; // LastParameter
        sup_bound[1] = self.curv.default_domain()[1];

        if self.first {
            inf_bound[2] = self.surf2.default_domain()[0];
            inf_bound[3] = self.surf2.default_domain()[2];
            sup_bound[2] = self.surf2.default_domain()[1];
            sup_bound[3] = self.surf2.default_domain()[3];
            if !is_infinite_value(inf_bound[2]) && !is_infinite_value(sup_bound[2]) {
                let range = sup_bound[2] - inf_bound[2];
                inf_bound[2] -= range;
                sup_bound[2] += range;
            }
            if !is_infinite_value(inf_bound[3]) && !is_infinite_value(sup_bound[3]) {
                let range = sup_bound[3] - inf_bound[3];
                inf_bound[3] -= range;
                sup_bound[3] += range;
            }
        } else {
            inf_bound[2] = self.surf1.default_domain()[0];
            inf_bound[3] = self.surf1.default_domain()[2];
            sup_bound[2] = self.surf1.default_domain()[1];
            sup_bound[3] = self.surf1.default_domain()[3];
            if !is_infinite_value(inf_bound[2]) && !is_infinite_value(sup_bound[2]) {
                let range = sup_bound[2] - inf_bound[2];
                inf_bound[2] -= range;
                sup_bound[2] += range;
            }
            if !is_infinite_value(inf_bound[3]) && !is_infinite_value(sup_bound[3]) {
                let range = sup_bound[3] - inf_bound[3];
                inf_bound[3] -= range;
                sup_bound[3] += range;
            }
        }
    }

    /// OCCT IsSolution(Sol, Tol) (BlendFunc_ChAsymInv.cxx L135-172).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        // OCCT: curv->D1(Sol(2), ptgui, d1gui); nplan = d1gui.Normalized();
        let ptgui = self.curv.point_at(sol[1]);
        let d1gui = self.curv.derivative_at(sol[1]);
        let nplan = d1gui.normalize_or_zero();
        let _ = ptgui;

        // OCCT: gp_Pnt2d pt2d(csurf->Value(Sol(1)));
        let csurf = self.csurf.expect("csurf");
        let pt2d = csurf.point_at(sol[0]);

        let (pts1, d1u1, d1v1, pts2);
        if self.first {
            // OCCT: surf1->D1(pt2d.X(), pt2d.Y(), pts1, d1u1, d1v1);
            let (p, du, dv) = self.surf1.derivatives(pt2d.x, pt2d.y);
            pts1 = p;
            d1u1 = du;
            d1v1 = dv;
            pts2 = self.surf2.point_at(sol[2], sol[3]);
        } else {
            // OCCT: surf1->D1(Sol(3), Sol(4), pts1, d1u1, d1v1);
            let (p, du, dv) = self.surf1.derivatives(sol[2], sol[3]);
            pts1 = p;
            d1u1 = du;
            d1v1 = dv;
            pts2 = self.surf2.point_at(pt2d.x, pt2d.y);
        }

        let nsurf1 = d1u1.cross(d1v1);
        let tsurf1 = nsurf1.cross(nplan);

        let s1s2 = pts2 - pts1;
        let psca_inv = 1.0 / tsurf1.dot(s1s2);
        let nordu1 = d1u1.length();
        let nordv1 = d1v1.length();

        let temp = 2.0 * (nordu1 + nordv1) * s1s2.length() + 2.0 * nordu1 * nordv1;

        let mut valsol = [0.0f64; 4];
        self.value(sol, &mut valsol);

        valsol[0].abs() < tol
            && valsol[1].abs() < tol
            && valsol[2].abs() < 2.0 * self.dist1 * tol
            && valsol[3].abs() < tol * (1.0 + self.tgang) * psca_inv.abs() * temp
    }

    /// OCCT ComputeValues(X, DegF, DegL)
    /// (BlendFunc_ChAsymInv.cxx L176-345).
    pub fn compute_values(&mut self, x: &[f64], deg_f: i32, deg_l: i32) -> bool {
        if deg_f > deg_l {
            return false;
        }

        let mut normg = 0.0f64;
        let pt2d: DVec2;
        let mut v2d = DVec2::ZERO;
        let nplan: DVec3;
        let mut dnplan = DVec3::ZERO;
        let ptgui: DVec3;
        let d1gui: DVec3;
        let pts1: DVec3;
        let pts2: DVec3;
        let d1u1: DVec3;
        let d1v1: DVec3;
        let d2u1: DVec3;
        let d2v1: DVec3;
        let d2uv1: DVec3;
        let d1u2: DVec3;
        let d1v2: DVec3;

        if (deg_f == 0) && (deg_l == 0) {
            // OCCT: curv->D1(X(2), ptgui, d1gui);
            let pg = self.curv.point_at(x[1]);
            let dg = self.curv.derivative_at(x[1]);
            nplan = dg.normalize_or_zero();
            ptgui = pg;
            d1gui = dg;

            // OCCT: pt2d = csurf->Value(X(1));
            let csurf = self.csurf.expect("csurf");
            pt2d = csurf.point_at(x[0]);

            if self.first {
                // OCCT: surf1->D1(pt2d.X(), pt2d.Y(), pts1, d1u1, d1v1);
                let (p, du, dv) = self.surf1.derivatives(pt2d.x, pt2d.y);
                pts1 = p;
                d1u1 = du;
                d1v1 = dv;
                pts2 = self.surf2.point_at(x[2], x[3]);
            } else {
                // OCCT: surf1->D1(X(3), X(4), pts1, d1u1, d1v1);
                let (p, du, dv) = self.surf1.derivatives(x[2], x[3]);
                pts1 = p;
                d1u1 = du;
                d1v1 = dv;
                pts2 = self.surf2.point_at(pt2d.x, pt2d.y);
            }
            d2u1 = DVec3::ZERO;
            d2v1 = DVec3::ZERO;
            d2uv1 = DVec3::ZERO;
            d1u2 = DVec3::ZERO;
            d1v2 = DVec3::ZERO;
        } else {
            // OCCT: curv->D2(X(2), ptgui, d1gui, d2gui);
            let pg = self.curv.point_at(x[1]);
            let dg = self.curv.derivative_at(x[1]);
            let d2g = self.curv.derivative2_at(x[1]);
            nplan = dg.normalize_or_zero();
            normg = dg.length();
            dnplan = (d2g - nplan.dot(d2g) * nplan) / normg;
            ptgui = pg;
            d1gui = dg;

            let csurf = self.csurf.expect("csurf");
            pt2d = csurf.point_at(x[0]);
            v2d = csurf.derivative_at(x[0]);

            if self.first {
                // OCCT: surf1->D2(pt2d.X(), pt2d.Y(), pts1, d1u1, d1v1, d2u1, d2v1, d2uv1);
                let (p, du, dv, d2u, d2uv, d2v) = self.surf1.derivatives2(pt2d.x, pt2d.y);
                pts1 = p;
                d1u1 = du;
                d1v1 = dv;
                d2u1 = d2u;
                d2uv1 = d2uv;
                d2v1 = d2v;
                // OCCT: surf2->D1(X(3), X(4), pts2, d1u2, d1v2);
                let (p2, du2, dv2) = self.surf2.derivatives(x[2], x[3]);
                pts2 = p2;
                d1u2 = du2;
                d1v2 = dv2;
            } else {
                // OCCT: surf1->D2(X(3), X(4), pts1, d1u1, d1v1, d2u1, d2v1, d2uv1);
                let (p, du, dv, d2u, d2uv, d2v) = self.surf1.derivatives2(x[2], x[3]);
                pts1 = p;
                d1u1 = du;
                d1v1 = dv;
                d2u1 = d2u;
                d2uv1 = d2uv;
                d2v1 = d2v;
                // OCCT: surf2->D1(pt2d.X(), pt2d.Y(), pts2, d1u2, d1v2);
                let (p2, du2, dv2) = self.surf2.derivatives(pt2d.x, pt2d.y);
                pts2 = p2;
                d1u2 = du2;
                d1v2 = dv2;
            }
        }

        // OCCT: gp_Vec nps1(ptgui, pts1), s1s2(pts1, pts2);
        let nps1 = pts1 - ptgui;
        let s1s2 = pts2 - pts1;
        let nsurf1 = d1u1.cross(d1v1);
        let tsurf1 = nsurf1.cross(nplan);
        let psca_inv = 1.0 / s1s2.dot(tsurf1);
        let f4 = nplan.dot(tsurf1.cross(s1s2)) * psca_inv;

        if deg_f == 0 {
            // OCCT: Dist = ptgui.XYZ().Dot(nplan.XYZ());
            let dist = ptgui.dot(nplan);
            self.fx[0] = pts1.dot(nplan) - dist;
            self.fx[1] = pts2.dot(nplan) - dist;
            self.fx[2] = self.dist1 * self.dist1 - nps1.length_squared();
            self.fx[3] = self.tgang - f4;
        }

        if deg_l == 1 {
            let nps2 = pts2 - ptgui;

            if self.first {
                // OCCT: dw1pts1 = v2d.X() * d1u1 + v2d.Y() * d1v1;
                let dw1pts1 = v2d.x * d1u1 + v2d.y * d1v1;
                let dw1du1 = v2d.x * d2u1 + v2d.y * d2uv1;
                let dw1dv1 = v2d.x * d2uv1 + v2d.y * d2v1;
                let dw1csurf = (dw1du1.cross(d1v1) + d1u1.cross(dw1dv1)).cross(nplan);
                let dwtsurf1 = nsurf1.cross(dnplan);

                self.dx[0][0] = nplan.dot(dw1pts1);
                self.dx[0][1] = dnplan.dot(nps1) - normg;
                self.dx[0][2] = 0.0;
                self.dx[0][3] = 0.0;

                self.dx[1][0] = 0.0;
                self.dx[1][1] = dnplan.dot(nps2) - normg;
                self.dx[1][2] = nplan.dot(d1u2);
                self.dx[1][3] = nplan.dot(d1v2);

                let temp_vec = 2.0 * nps1;
                self.dx[2][0] = -dw1pts1.dot(temp_vec);
                self.dx[2][1] = d1gui.dot(temp_vec);
                self.dx[2][2] = 0.0;
                self.dx[2][3] = 0.0;

                let mut temp = f4 * (dw1csurf.dot(s1s2) - tsurf1.dot(dw1pts1));
                temp += nplan.dot(tsurf1.cross(dw1pts1) - dw1csurf.cross(s1s2));
                self.dx[3][0] = psca_inv * temp;

                let mut temp = f4 * dwtsurf1.dot(s1s2);
                temp -= dnplan.dot(temp_vec) + nplan.dot(dwtsurf1.cross(s1s2));
                self.dx[3][1] = psca_inv * temp;
                let temp = f4 * tsurf1.dot(d1u2) - nplan.dot(tsurf1.cross(d1u2));
                self.dx[3][2] = psca_inv * temp;

                let temp = f4 * tsurf1.dot(d1v2) - nplan.dot(tsurf1.cross(d1v2));
                self.dx[3][3] = psca_inv * temp;
            } else {
                // OCCT: d1utsurf1 = (d2u1.Crossed(d1v1) + d1u1.Crossed(d2uv1)).Crossed(nplan);
                let d1utsurf1 = (d2u1.cross(d1v1) + d1u1.cross(d2uv1)).cross(nplan);
                let d1vtsurf1 = (d2uv1.cross(d1v1) + d1u1.cross(d2v1)).cross(nplan);
                let dw2pts2 = v2d.x * d1u2 + v2d.y * d1v2;
                let dwtsurf1 = nsurf1.cross(dnplan);

                self.dx[0][0] = 0.0;
                self.dx[0][1] = dnplan.dot(nps1) - normg;
                self.dx[0][2] = nplan.dot(d1u1);
                self.dx[0][3] = nplan.dot(d1v1);

                self.dx[1][0] = nplan.dot(dw2pts2);
                self.dx[1][1] = dnplan.dot(nps2) - normg;
                self.dx[1][2] = 0.0;
                self.dx[1][3] = 0.0;

                let mut temp_vec = 2.0 * nps1;
                self.dx[2][0] = 0.0;
                self.dx[2][1] = d1gui.dot(temp_vec);

                temp_vec = -temp_vec;
                self.dx[2][2] = d1u1.dot(temp_vec);
                self.dx[2][3] = d1v1.dot(temp_vec);

                let temp = f4 * tsurf1.dot(dw2pts2) - nplan.dot(tsurf1.cross(dw2pts2));
                self.dx[3][0] = psca_inv * temp;

                let mut temp = f4 * dwtsurf1.dot(s1s2);
                temp -= dnplan.dot(temp_vec) + nplan.dot(dwtsurf1.cross(s1s2));
                self.dx[3][1] = psca_inv * temp;

                let mut temp = f4 * (d1utsurf1.dot(s1s2) - tsurf1.dot(d1u1));
                temp += nplan.dot(tsurf1.cross(d1u1) - d1utsurf1.cross(s1s2));
                self.dx[3][2] = psca_inv * temp;

                let mut temp = f4 * (d1vtsurf1.dot(s1s2) - tsurf1.dot(d1v1));
                temp += nplan.dot(tsurf1.cross(d1v1) - d1vtsurf1.cross(s1s2));
                self.dx[3][3] = psca_inv * temp;
            }
        }

        true
    }

    /// OCCT Value(X, F) (BlendFunc_ChAsymInv.cxx L349-354).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let error = self.compute_values(x, 0, 0);
        f.copy_from_slice(&self.fx);
        error
    }

    /// OCCT Derivatives(X, D) (BlendFunc_ChAsymInv.cxx L358-363).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        let error = self.compute_values(x, 1, 1);
        // OCCT: D = DX;
        for (drow, xrow) in d.iter_mut().zip(self.dx.iter()) {
            drow.clone_from(xrow);
        }
        error
    }

    /// OCCT Values(X, F, D) (BlendFunc_ChAsymInv.cxx L367-373).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        let error = self.compute_values(x, 0, 1);
        f.copy_from_slice(&self.fx);
        // OCCT: D = DX;
        for (drow, xrow) in d.iter_mut().zip(self.dx.iter()) {
            drow.clone_from(xrow);
        }
        error
    }
}

impl<'a> FunctionSetWithDerivatives for BlendFuncChAsymInv<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: inherited from Blend_FuncInv::NbVariables (returns 4).
        BlendFuncInv::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendFuncChAsymInv::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendFuncChAsymInv::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncChAsymInv::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncChAsymInv::values(self, x, f, df)
    }
}

impl<'a> BlendFuncInv for BlendFuncChAsymInv<'a> {
    fn set_curve_on_surface(&mut self, on_first: bool, c_on_surf: &Curve2d) {
        // OCCT ChAsymInv.cxx L60-64; the rcad port stores the curve reference
        // (OCCT copies the handle).  SAFETY: the caller owns the Curve2d for
        // the lifetime 'a of this function object — same invariant as the
        // OCCT handle (the restriction curve outlives the solver).
        let c: &'a Curve2d = unsafe { &*(c_on_surf as *const Curve2d) };
        BlendFuncChAsymInv::set_curve_on_surface(self, on_first, c)
    }

    fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        BlendFuncChAsymInv::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendFuncChAsymInv::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendFuncChAsymInv::is_solution(self, sol, tol)
    }
}
