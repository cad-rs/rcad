//! OCCT BlendFunc_EvolRad (TKFillet/BlendFunc) — 1:1 port of
//! BlendFunc_EvolRad.hxx (L38-282) + BlendFunc_EvolRad.cxx (whole file
//! L39-2036).  `BRepBlend_EvolRad` is a pure typedef
//! (`typedef BlendFunc_EvolRad BRepBlend_EvolRad;`, BRepBlend_EvolRad.hxx
//! L21) so this single translation serves both names.
//!
//! [`BlendFuncEvolRad`]: the constructor, the Set overloads, IsSolution, the
//! accessors, the interval family and the three solver/approx trait impls
//! live here; ComputeValues and the four Section overloads live in
//! [`super::brep_blend_func_evolrad_b`].
//!
//! Architecture mappings (mirroring [`super::brep_blend_func_consrad`]):
//! `class BlendFunc_EvolRad : public Blend_Function` is expressed by
//! implementing the [`BlendFunction`] trait over the
//! `math_FunctionSetWithDerivatives` base; `tcurv` (OCCT handle aliasing
//! `curv` until Set(First, Last)) maps to an owned trimmed [`Curve3`] copy
//! plus an accessor, and likewise `tevol` (aliasing `fevol` until Set(First,
//! Last) trims it) maps to `Option<LawFunctionHandle>`; `math_Vector`/
//! `math_Matrix` map to `[f64; 4]` / `Vec<Vec<f64>>`; `BlendFunc_Tensor`
//! is reused from [`super::brep_blend_func_consrad::BlendFuncTensor`];
//! `gp_Circ` maps to the kernel [`Circle3`].
//! Pending kernel dependency (marked GAP, plan 0.6): math_SVD (the
//! second-chance solver of the D1/D2 Sections — the OCCT `!IsDone()` path
//! is preserved).

use rcad_kernel::core::precision::is_infinite_value;
use rcad_kernel::geom::{Curve3, Surface3, SurfaceEval as _};
use rcad_kernel::math::function_set_root::FunctionSetWithDerivatives;
use rcad_kernel::math::{GeomAbsShape, MatD, VecD};

use crate::geomalgo::law::LawFunctionHandle;

use super::brep_blend_func::{
    blend_func_get_minimal_weights, blend_func_get_shape, blend_func_next_shape,
    BlendFuncSectionShape, ConvertParameterisationType,
};
use super::brep_blend_function::{BlendAppFunction, BlendFunction};
use super::brep_blend_func_consrad::{
    geomfill_get_tolerance, geomfill_knots, geomfill_mults, BlendFuncTensor,
};

/// OCCT Eps constant (BlendFunc_EvolRad.cxx L37).
pub(crate) const EPS: f64 = 1.0e-15;

/// OCCT BlendFunc_EvolRad — function for a variable-radius rolling-ball
/// blending surface; the radius evolution law is a Law_Function
/// (BlendFunc_EvolRad.hxx L38).
pub struct BlendFuncEvolRad<'a> {
    // OCCT BlendFunc_EvolRad.hxx fields (L206-281).
    pub(crate) surf1: &'a Surface3,
    pub(crate) surf2: &'a Surface3,
    pub(crate) curv: &'a Curve3,
    // OCCT: occ::handle<Adaptor3d_Curve> tcurv — a handle aliasing curv until
    // Set(First, Last) trims it.  The rcad port owns the trimmed copy.
    pub(crate) tcurv: Option<Curve3>,
    // OCCT: occ::handle<Law_Function> fevol — the radius evolution law.
    pub(crate) fevol: LawFunctionHandle,
    // OCCT: occ::handle<Law_Function> tevol — a handle aliasing fevol until
    // Set(First, Last) trims it.
    pub(crate) tevol: Option<LawFunctionHandle>,
    pub(crate) pts1: DVec3,
    pub(crate) pts2: DVec3,
    pub(crate) istangent: bool,
    pub(crate) tg1: DVec3,
    pub(crate) tg12d: DVec2,
    pub(crate) tg2: DVec3,
    pub(crate) tg22d: DVec2,
    pub(crate) param: f64,
    pub(crate) sg1: f64,
    pub(crate) sg2: f64,
    pub(crate) ray: f64,
    pub(crate) dray: f64,
    pub(crate) d2ray: f64,
    pub(crate) choix: i32,
    pub(crate) my_x_order: i32,
    pub(crate) my_t_order: i32,
    pub(crate) xval: [f64; 4],
    pub(crate) tval: f64,
    pub(crate) d1u1: DVec3,
    pub(crate) d1u2: DVec3,
    pub(crate) d1v1: DVec3,
    pub(crate) d1v2: DVec3,
    pub(crate) d2u1: DVec3,
    pub(crate) d2v1: DVec3,
    pub(crate) d2uv1: DVec3,
    pub(crate) d2u2: DVec3,
    pub(crate) d2v2: DVec3,
    pub(crate) d2uv2: DVec3,
    pub(crate) dn1w: DVec3,
    pub(crate) dn2w: DVec3,
    pub(crate) d2n1w: DVec3,
    pub(crate) d2n2w: DVec3,
    pub(crate) nplan: DVec3,
    pub(crate) nsurf1: DVec3,
    pub(crate) nsurf2: DVec3,
    pub(crate) dns1u1: DVec3,
    pub(crate) dns1u2: DVec3,
    pub(crate) dns1v1: DVec3,
    pub(crate) dns1v2: DVec3,
    pub(crate) dnplan: DVec3,
    pub(crate) d2nplan: DVec3,
    // OCCT BlendFunc_EvolRad.hxx L253-254 — declared but never accessed by
    // the EvolRad bodies (vestigial fields, kept for structure parity).
    #[allow(dead_code)]
    pub(crate) dnsurf1: DVec3,
    #[allow(dead_code)]
    pub(crate) dnsurf2: DVec3,
    pub(crate) dndu1: DVec3,
    pub(crate) dndu2: DVec3,
    pub(crate) dndv1: DVec3,
    pub(crate) dndv2: DVec3,
    pub(crate) d2ndu1: DVec3,
    pub(crate) d2ndu2: DVec3,
    pub(crate) d2ndv1: DVec3,
    pub(crate) d2ndv2: DVec3,
    pub(crate) d2nduv1: DVec3,
    pub(crate) d2nduv2: DVec3,
    pub(crate) d2ndtu1: DVec3,
    pub(crate) d2ndtu2: DVec3,
    pub(crate) d2ndtv1: DVec3,
    pub(crate) d2ndtv2: DVec3,
    // OCCT ComputeValues `static` locals (EvolRad.cxx L201-204) that must
    // survive across calls (read on cached t-states): the guide point /
    // derivatives, 1/normtg and d(1/normtg)/dt.
    pub(crate) ptgui: DVec3,
    pub(crate) d1gui: DVec3,
    pub(crate) d2gui: DVec3,
    pub(crate) invnormtg: f64,
    pub(crate) dinvnormtg: f64,
    pub(crate) e: [f64; 4],
    pub(crate) dedx: Vec<Vec<f64>>,
    pub(crate) dedt: [f64; 4],
    pub(crate) d2edx2: BlendFuncTensor,
    pub(crate) d2edxdt: Vec<Vec<f64>>,
    pub(crate) d2edt2: [f64; 4],
    pub(crate) minang: f64,
    pub(crate) maxang: f64,
    pub(crate) lengthmin: f64,
    pub(crate) lengthmax: f64,
    pub(crate) distmin: f64,
    pub(crate) my_s_shape: BlendFuncSectionShape,
    pub(crate) my_t_conv: ConvertParameterisationType,
}

use glam::{DVec2, DVec3};

impl<'a> BlendFuncEvolRad<'a> {
    /// Architecture mapping: OCCT `tcurv` is a handle aliasing `curv` until
    /// Set(First, Last) replaces it with a trimmed copy (EvolRad.cxx L860).
    #[inline]
    pub(crate) fn tcurv(&self) -> &Curve3 {
        self.tcurv.as_ref().unwrap_or(self.curv)
    }

    /// Architecture mapping: OCCT `tevol` is a handle aliasing `fevol` until
    /// Set(First, Last) replaces it with a trimmed law (EvolRad.cxx L861).
    #[inline]
    pub(crate) fn tevol(&self) -> &LawFunctionHandle {
        self.tevol.as_ref().unwrap_or(&self.fevol)
    }

    /// OCCT BlendFunc_EvolRad(S1, S2, C, Law) (BlendFunc_EvolRad.cxx
    /// L102-133).
    pub fn new(
        s1: &'a Surface3,
        s2: &'a Surface3,
        c: &'a Curve3,
        law: LawFunctionHandle,
    ) -> Self {
        let mut func = BlendFuncEvolRad {
            surf1: s1,
            surf2: s2,
            curv: c,
            tcurv: None, // OCCT: tcurv(C) — alias, see tcurv().
            fevol: law,
            tevol: None, // OCCT: tevol(Law) — alias, see tevol().
            pts1: DVec3::ZERO,
            pts2: DVec3::ZERO,
            istangent: true,
            tg1: DVec3::ZERO,
            tg12d: DVec2::ZERO,
            tg2: DVec3::ZERO,
            tg22d: DVec2::ZERO,
            param: 0.0,
            sg1: 0.0,
            sg2: 0.0,
            ray: 0.0,
            dray: 0.0,
            d2ray: 0.0,
            choix: 0,
            my_x_order: -1,
            my_t_order: -1,
            xval: [-9.876e100; 4], // OCCT: xval(1, 4) then xval.Init(-9.876e100)
            tval: -9.876e100,      // OCCT: tval = -9.876e100
            d1u1: DVec3::ZERO,
            d1u2: DVec3::ZERO,
            d1v1: DVec3::ZERO,
            d1v2: DVec3::ZERO,
            d2u1: DVec3::ZERO,
            d2v1: DVec3::ZERO,
            d2uv1: DVec3::ZERO,
            d2u2: DVec3::ZERO,
            d2v2: DVec3::ZERO,
            d2uv2: DVec3::ZERO,
            dn1w: DVec3::ZERO,
            dn2w: DVec3::ZERO,
            d2n1w: DVec3::ZERO,
            d2n2w: DVec3::ZERO,
            nplan: DVec3::ZERO,
            nsurf1: DVec3::ZERO,
            nsurf2: DVec3::ZERO,
            dns1u1: DVec3::ZERO,
            dns1u2: DVec3::ZERO,
            dns1v1: DVec3::ZERO,
            dns1v2: DVec3::ZERO,
            dnplan: DVec3::ZERO,
            d2nplan: DVec3::ZERO,
            dnsurf1: DVec3::ZERO,
            dnsurf2: DVec3::ZERO,
            dndu1: DVec3::ZERO,
            dndu2: DVec3::ZERO,
            dndv1: DVec3::ZERO,
            dndv2: DVec3::ZERO,
            d2ndu1: DVec3::ZERO,
            d2ndu2: DVec3::ZERO,
            d2ndv1: DVec3::ZERO,
            d2ndv2: DVec3::ZERO,
            d2nduv1: DVec3::ZERO,
            d2nduv2: DVec3::ZERO,
            d2ndtu1: DVec3::ZERO,
            d2ndtu2: DVec3::ZERO,
            d2ndtv1: DVec3::ZERO,
            d2ndtv2: DVec3::ZERO,
            ptgui: DVec3::ZERO,
            d1gui: DVec3::ZERO,
            d2gui: DVec3::ZERO,
            invnormtg: 0.0,
            dinvnormtg: 0.0,
            e: [0.0; 4],
            dedx: vec![vec![0.0; 4]; 4],
            dedt: [0.0; 4],
            d2edx2: BlendFuncTensor::new(4, 4, 4),
            d2edxdt: vec![vec![0.0; 4]; 4],
            d2edt2: [0.0; 4],
            minang: f64::MAX,  // OCCT: RealLast()
            maxang: -f64::MAX, // OCCT: RealFirst()
            lengthmin: f64::MAX,   // OCCT: RealLast()
            lengthmax: -f64::MAX,  // OCCT: RealFirst()
            distmin: f64::MAX, // OCCT: RealLast()
            my_s_shape: BlendFuncSectionShape::Rational,
            my_t_conv: ConvertParameterisationType::TgtThetaOver2,
        };
        // OCCT L125-126: fevol = Law; tevol = Law (handled in the initializer
        // through the shared handle: tevol aliases fevol, see tevol()).
        // OCCT L129-132: initialisation of the cash control variables.
        func.tval = -9.876e100;
        func.xval = [-9.876e100; 4];
        func.my_x_order = -1;
        func.my_t_order = -1;
        func
    }

    /// OCCT NbEquations() (BlendFunc_EvolRad.cxx L137-140) — returns 4.
    pub fn nb_equations(&self) -> usize {
        4
    }

    /// OCCT Set(Choix) (BlendFunc_EvolRad.cxx L144-176) — inits the
    /// "quadrant" signs.
    pub fn set(&mut self, choix: i32) {
        self.choix = choix;
        match self.choix {
            1 | 2 => {
                self.sg1 = -1.0;
                self.sg2 = -1.0;
            }
            3 | 4 => {
                self.sg1 = 1.0;
                self.sg2 = -1.0;
            }
            5 | 6 => {
                self.sg1 = 1.0;
                self.sg2 = 1.0;
            }
            7 | 8 => {
                self.sg1 = -1.0;
                self.sg2 = 1.0;
            }
            _ => {
                self.sg1 = -1.0;
                self.sg2 = -1.0;
            }
        }
    }

    /// OCCT Set(TypeSection) (BlendFunc_EvolRad.cxx L180-183) — sets the
    /// type of section generation for the approximations.
    pub fn set_section_shape(&mut self, type_section: BlendFuncSectionShape) {
        self.my_s_shape = type_section;
    }

    /// OCCT Set(Param) (BlendFunc_EvolRad.cxx L847-850).
    pub fn set_param(&mut self, param: f64) {
        self.param = param;
    }

    /// OCCT Set(First, Last) (BlendFunc_EvolRad.cxx L858-862) — segments the
    /// curve and the law in their useful part.
    pub fn set_interval(&mut self, first: f64, last: f64) {
        // OCCT: tcurv = curv->Trim(First, Last, 1.e-12);
        self.tcurv = Some(Curve3::Trimmed(rcad_kernel::geom::TrimmedCurve3 {
            curve: Box::new(self.curv.clone()),
            first,
            last,
        }));
        // OCCT: tevol = fevol->Trim(First, Last, 1.e-12);
        self.tevol = Some(self.fevol.borrow().trim(first, last, 1.0e-12));
    }

    /// OCCT GetTolerance(Tolerance, Tol) (BlendFunc_EvolRad.cxx L866-872).
    pub fn get_tolerance(&self, tolerance: &mut [f64], tol: f64) {
        tolerance[0] = self.surf1.u_resolution(tol);
        tolerance[1] = self.surf1.v_resolution(tol);
        tolerance[2] = self.surf2.u_resolution(tol);
        tolerance[3] = self.surf2.v_resolution(tol);
    }

    /// OCCT GetBounds(InfBound, SupBound) (BlendFunc_EvolRad.cxx L876-896).
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

    /// OCCT IsSolution(Sol, Tol) (BlendFunc_EvolRad.cxx L900-1018).
    pub fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        let ok = self.compute_values(sol, 1, true, self.param);

        if self.e[0].abs() <= tol
            && self.e[1] * self.e[1] + self.e[2] * self.e[2] + self.e[3] * self.e[3] <= tol * tol
        {
            // ns1, ns2, np are copied locally to avoid crushing the fields !
            let mut ns1 = self.nsurf1;
            let mut ns2 = self.nsurf2;
            let np = self.nplan;

            let mut norm = self.nplan.cross(ns1).length();
            if norm < EPS {
                norm = 1.0; // Unsatisfactory, but it is not necessary to stop
            }
            // OCCT: ns1.SetLinearForm(nplan.Dot(ns1) / norm, nplan, -1. / norm, ns1);
            ns1 = (self.nplan.dot(ns1) / norm) * self.nplan + (-1.0 / norm) * ns1;

            let mut norm = self.nplan.cross(ns2).length();
            if norm < EPS {
                norm = 1.0; // Unsatisfactory, but it is not necessary to stop
            }
            ns2 = (self.nplan.dot(ns2) / norm) * self.nplan + (-1.0 / norm) * ns2;

            // OCCT: double maxpiv = 1.e-14; math_Gauss Resol(DEDX, maxpiv).
            // (The rcad MathGauss models the pivoted Gauss solve; the pivot
            // floor parameter has no rcad equivalent.)
            let mut a = MatD::new(4, 4);
            for r in 1..=4 {
                for c in 1..=4 {
                    a.set(r, c, self.dedx[r - 1][c - 1]);
                }
            }
            self.istangent = false;
            // OCCT: math_Vector controle(1, 4), solution(1, 4), tolerances(1, 4);
            //       GetTolerance(tolerances, Tol);
            let mut tolerances = [0.0f64; 4];
            self.get_tolerance(&mut tolerances, tol);
            // OCCT: math_Vector solution(1, 4) — the solver output shared by
            // the Gauss branch.
            let mut solution = [0.0f64; 4];
            let resol = rcad_kernel::math::math_gauss::MathGauss::new(&a);
            if resol.is_done() {
                // OCCT: Resol.Solve(-DEDT, solution);
                let mut x = VecD::new(4);
                for i in 1..=4 {
                    x.set(i, -self.dedt[i - 1]);
                }
                resol.solve(&mut x);
                for i in 1..=4 {
                    solution[i - 1] = x.get(i);
                }
                // OCCT: controle = DEDT.Added(DEDX.Multiplied(solution));
                let mut controle = [0.0f64; 4];
                for i in 1..=4 {
                    let mut somme = 0.0;
                    for j in 1..=4 {
                        somme += self.dedx[i - 1][j - 1] * solution[j - 1];
                    }
                    controle[i - 1] = self.dedt[i - 1] + somme;
                }
                if controle[0].abs() > tolerances[0]
                    || controle[1].abs() > tolerances[1]
                    || controle[2].abs() > tolerances[2]
                    || controle[3].abs() > tolerances[3]
                {
                    self.istangent = true;
                }
                if !self.istangent {
                    // OCCT: tg1.SetLinearForm(solution(1), d1u1, solution(2), d1v1);
                    self.tg1 = solution[0] * self.d1u1 + solution[1] * self.d1v1;
                    self.tg2 = solution[2] * self.d1u2 + solution[3] * self.d1v2;
                    self.tg12d = DVec2::new(solution[0], solution[1]);
                    self.tg22d = DVec2::new(solution[2], solution[3]);
                }
            } else {
                self.istangent = true;
            }

            // update of maxang

            if self.sg1 > 0.0 {
                // sg1*ray
                ns1 = -ns1;
            }
            if self.sg2 > 0.0 {
                // sg2*ray
                ns2 = -ns2;
            }
            let mut cosa = ns1.dot(ns2);
            let mut sina = np.dot(ns1.cross(ns2));
            if self.choix % 2 != 0 {
                sina = -sina; // nplan is changed into -nplan
            }

            if cosa > 1.0 {
                cosa = 1.0;
                sina = 0.0;
            }
            let mut angle = cosa.acos();
            // Reframing on ]-pi/2, 3pi/2]
            if sina < 0.0 {
                if cosa > 0.0 {
                    angle = -angle;
                } else {
                    angle = 2.0 * std::f64::consts::PI - angle;
                }
            }

            if angle.abs() > self.maxang {
                self.maxang = angle.abs();
            }
            if angle.abs() < self.minang {
                self.minang = angle.abs();
            }
            if (angle * self.ray).abs() < self.lengthmin {
                self.lengthmin = (angle * self.ray).abs();
            }
            if (angle * self.ray).abs() > self.lengthmax {
                self.lengthmax = (angle * self.ray).abs();
            }
            self.distmin = self.distmin.min(self.pts1.distance(self.pts2));

            return ok;
        }
        self.istangent = true;
        false
    }

    /// OCCT GetMinimalDistance() (BlendFunc_EvolRad.cxx L1022-1025).
    pub fn get_minimal_distance(&self) -> f64 {
        self.distmin
    }

    /// OCCT Value(X, F) (BlendFunc_EvolRad.cxx L1029-1035).
    pub fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        let ok = self.compute_values(x, 0, false, 0.0);
        f.copy_from_slice(&self.e);
        ok
    }

    /// OCCT Derivatives(X, D) (BlendFunc_EvolRad.cxx L1039-1045).
    pub fn derivatives(&mut self, x: &[f64], d: &mut [Vec<f64>]) -> bool {
        let ok = self.compute_values(x, 1, false, 0.0);
        for i in 0..4 {
            d[i].copy_from_slice(&self.dedx[i]);
        }
        ok
    }

    /// OCCT Values(X, F, D) (BlendFunc_EvolRad.cxx L1047-1054).
    pub fn values(&mut self, x: &[f64], f: &mut [f64], d: &mut [Vec<f64>]) -> bool {
        let ok = self.compute_values(x, 1, false, 0.0);
        f.copy_from_slice(&self.e);
        for i in 0..4 {
            d[i].copy_from_slice(&self.dedx[i]);
        }
        ok
    }

    /// OCCT PointOnS1() (BlendFunc_EvolRad.cxx L1197-1200).
    pub fn point_on_s1(&self) -> DVec3 {
        self.pts1
    }

    /// OCCT PointOnS2() (BlendFunc_EvolRad.cxx L1204-1207).
    pub fn point_on_s2(&self) -> DVec3 {
        self.pts2
    }

    /// OCCT IsTangencyPoint() (BlendFunc_EvolRad.cxx L1211-1214).
    pub fn is_tangency_point(&self) -> bool {
        self.istangent
    }

    /// OCCT TangentOnS1() (BlendFunc_EvolRad.cxx L1218-1225).
    pub fn tangent_on_s1(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError");
        }
        self.tg1
    }

    /// OCCT TangentOnS2() (BlendFunc_EvolRad.cxx L1229-1236).
    pub fn tangent_on_s2(&self) -> DVec3 {
        if self.istangent {
            panic!("Standard_DomainError");
        }
        self.tg2
    }

    /// OCCT Tangent2dOnS1() (BlendFunc_EvolRad.cxx L1240-1247).
    pub fn tangent_2d_on_s1(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError");
        }
        self.tg12d
    }

    /// OCCT Tangent2dOnS2() (BlendFunc_EvolRad.cxx L1251-1258).
    pub fn tangent_2d_on_s2(&self) -> DVec2 {
        if self.istangent {
            panic!("Standard_DomainError");
        }
        self.tg22d
    }

    /// OCCT TwistOnS1() (BlendFunc_EvolRad.cxx L1113-1120).
    pub fn twist_on_s1(&self) -> bool {
        if self.istangent {
            panic!("Standard_DomainError");
        }
        self.tg1.dot(self.nplan) < 0.0
    }

    /// OCCT TwistOnS2() (BlendFunc_EvolRad.cxx L1124-1131).
    pub fn twist_on_s2(&self) -> bool {
        if self.istangent {
            panic!("Standard_DomainError");
        }
        self.tg2.dot(self.nplan) < 0.0
    }

    /// OCCT IsRational() (BlendFunc_EvolRad.cxx L1262-1265).
    pub fn is_rational(&self) -> bool {
        self.my_s_shape == BlendFuncSectionShape::Rational
            || self.my_s_shape == BlendFuncSectionShape::QuasiAngular
    }

    /// OCCT GetSectionSize() (BlendFunc_EvolRad.cxx L1269-1272).
    pub fn get_section_size(&self) -> f64 {
        self.lengthmax
    }

    /// OCCT GetMinimalWeight(Weights) (BlendFunc_EvolRad.cxx L1276-1279).
    pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
        blend_func_get_minimal_weights(self.my_s_shape, self.my_t_conv, self.minang, self.maxang, weigths);
    }

    /// OCCT NbIntervals(S) (BlendFunc_EvolRad.cxx L1283-1302).
    pub fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        let nb_int_courbe = self.curv.nb_intervals(blend_func_next_shape(s));
        let nb_int_loi = self.fevol.borrow().nb_intervals(s);

        if nb_int_loi == 1 {
            return nb_int_courbe;
        }

        let mut int_c = vec![0.0f64; nb_int_courbe + 1];
        let mut int_l = vec![0.0f64; nb_int_loi + 1];
        let mut inter: Vec<f64> = Vec::new();
        self.curv.intervals(&mut int_c, blend_func_next_shape(s));
        self.fevol.borrow().intervals(&mut int_l, s);

        fusionne_intervalles(&int_c, &int_l, &mut inter);
        inter.len() - 1
    }

    /// OCCT Intervals(T, S) (BlendFunc_EvolRad.cxx L1306-1330).
    pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let nb_int_courbe = self.curv.nb_intervals(blend_func_next_shape(s));
        let nb_int_loi = self.fevol.borrow().nb_intervals(s);

        if nb_int_loi == 1 {
            let mut intervals = Vec::new();
            self.curv.intervals(&mut intervals, blend_func_next_shape(s));
            for (dst, src) in t.iter_mut().zip(intervals) {
                *dst = src;
            }
        } else {
            let mut int_c = vec![0.0f64; nb_int_courbe + 1];
            let mut int_l = vec![0.0f64; nb_int_loi + 1];
            let mut inter: Vec<f64> = Vec::new();
            self.curv.intervals(&mut int_c, blend_func_next_shape(s));
            self.fevol.borrow().intervals(&mut int_l, s);

            fusionne_intervalles(&int_c, &int_l, &mut inter);
            for (dst, src) in t.iter_mut().zip(inter) {
                *dst = src;
            }
        }
    }

    /// OCCT GetShape(NbPoles, NbKnots, Degree, NbPoles2d)
    /// (BlendFunc_EvolRad.cxx L1334-1338).
    pub fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        *nb_poles_2d = 2;
        blend_func_get_shape(
            self.my_s_shape,
            self.maxang,
            nb_poles,
            nb_knots,
            degree,
            &mut self.my_t_conv,
        );
    }

    /// OCCT GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d, Tol1d)
    /// (BlendFunc_EvolRad.cxx L1344-1358) — tolerances used for
    /// approximations.
    pub fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        let low = 0usize; // OCCT: Tol3d.Lower()
        let up = tol3d.len() - 1; // OCCT: Tol3d.Upper()
        let rayon = self.lengthmin / self.maxang; // a radius is subtracted
        let tol = geomfill_get_tolerance(self.my_t_conv, self.maxang, rayon, angle_tol, surf_tol);
        for v in tol1d.iter_mut() {
            *v = surf_tol;
        }
        for v in tol3d.iter_mut() {
            *v = surf_tol;
        }
        tol3d[low + 1] = tol.min(surf_tol);
        tol3d[up - 1] = tol.min(surf_tol);
        tol3d[low] = tol.min(bound_tol);
        tol3d[up] = tol.min(bound_tol);
    }

    /// OCCT Knots(TKnots) (BlendFunc_EvolRad.cxx L1362-1365).
    pub fn knots(&mut self, tknots: &mut [f64]) {
        geomfill_knots(self.my_t_conv, tknots);
    }

    /// OCCT Mults(TMults) (BlendFunc_EvolRad.cxx L1369-1372).
    pub fn mults(&mut self, tmults: &mut [i32]) {
        geomfill_mults(self.my_t_conv, tmults);
    }

    /// OCCT Resolution(IC2d, Tol, TolU, TolV) (BlendFunc_EvolRad.cxx
    /// L2021-2036).
    pub fn resolution(&self, ic_2d: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        if ic_2d == 1 {
            *tol_u = self.surf1.u_resolution(tol);
            *tol_v = self.surf1.v_resolution(tol);
        } else {
            *tol_u = self.surf2.u_resolution(tol);
            *tol_v = self.surf2.v_resolution(tol);
        }
    }
}

/// OCCT FusionneIntervalles(I1, I2, Seq) (BlendFunc_EvolRad.cxx L39-98) —
/// merges the two interval tables, removing multiple occurrences.
pub(crate) fn fusionne_intervalles(i1: &[f64], i2: &[f64], seq: &mut Vec<f64>) {
    let mut ind1 = 0usize; // OCCT: int ind1 = 1 (1-based)
    let mut ind2 = 0usize; // OCCT: int ind2 = 1
    let epspar = 0.99 * rcad_kernel::core::precision::p_confusion();
    // supposed that positioning works with PConfusion()/2

    //--- TABSOR is filled by parsing TABLE1 and TABLE2 simultaneously ---
    //------------------ by removing multiple occurrences ------------

    while ind1 < i1.len() && ind2 < i2.len() {
        let v1 = i1[ind1];
        let v2 = i2[ind2];
        if (v1 - v2).abs() <= epspar {
            // Here elements of I1 and I2 are suitable.
            seq.push((v1 + v2) / 2.0);
            ind1 += 1;
            ind2 += 1;
        } else if v1 < v2 {
            // Here the element of I1 is suitable.
            seq.push(v1);
            ind1 += 1;
        } else {
            // Here the element of TABLE2 is suitable.
            seq.push(v2);
            ind2 += 1;
        }
    }

    if ind1 >= i1.len() {
        //----- Here I1 is empty, to be completed with the end of TABLE2 -------
        while ind2 < i2.len() {
            seq.push(i2[ind2]);
            ind2 += 1;
        }
    }

    if ind2 >= i2.len() {
        //----- Here I2 is empty, to be completed with the end of I1 -------
        while ind1 < i1.len() {
            seq.push(i1[ind1]);
            ind1 += 1;
        }
    }
}

impl<'a> FunctionSetWithDerivatives for BlendFuncEvolRad<'a> {
    fn nb_variables(&self) -> usize {
        // OCCT: inherited from Blend_Function::NbVariables (returns 4).
        BlendFunction::nb_variables(self)
    }

    fn nb_equations(&self) -> usize {
        BlendFuncEvolRad::nb_equations(self)
    }

    fn value(&mut self, x: &[f64], f: &mut [f64]) -> bool {
        BlendFuncEvolRad::value(self, x, f)
    }

    fn derivatives(&mut self, x: &[f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncEvolRad::derivatives(self, x, df)
    }

    fn values(&mut self, x: &[f64], f: &mut [f64], df: &mut [Vec<f64>]) -> bool {
        BlendFuncEvolRad::values(self, x, f, df)
    }
}

impl<'a> BlendAppFunction for BlendFuncEvolRad<'a> {
    fn set_param(&mut self, param: f64) {
        BlendFuncEvolRad::set_param(self, param)
    }

    fn set_interval(&mut self, first: f64, last: f64) {
        BlendFuncEvolRad::set_interval(self, first, last)
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
        BlendFuncEvolRad::get_tolerance(self, tolerance, tol)
    }

    fn get_bounds(&self, inf_bound: &mut [f64], sup_bound: &mut [f64]) {
        BlendFuncEvolRad::get_bounds(self, inf_bound, sup_bound)
    }

    fn is_solution(&mut self, sol: &[f64], tol: f64) -> bool {
        BlendFuncEvolRad::is_solution(self, sol, tol)
    }

    fn get_minimal_distance(&self) -> f64 {
        BlendFuncEvolRad::get_minimal_distance(self)
    }

    fn is_rational(&self) -> bool {
        BlendFuncEvolRad::is_rational(self)
    }

    fn get_section_size(&self) -> f64 {
        BlendFuncEvolRad::get_section_size(self)
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        BlendFuncEvolRad::get_minimal_weight(self, weigths)
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> usize {
        BlendFuncEvolRad::nb_intervals(self, s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        BlendFuncEvolRad::intervals(self, t, s)
    }

    fn get_shape(
        &mut self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles_2d: &mut i32,
    ) {
        BlendFuncEvolRad::get_shape(self, nb_poles, nb_knots, degree, nb_poles_2d)
    }

    fn get_approx_tolerance(
        &self,
        bound_tol: f64,
        surf_tol: f64,
        angle_tol: f64,
        tol3d: &mut [f64],
        tol1d: &mut [f64],
    ) {
        BlendFuncEvolRad::get_approx_tolerance(self, bound_tol, surf_tol, angle_tol, tol3d, tol1d)
    }

    fn knots(&mut self, tknots: &mut [f64]) {
        BlendFuncEvolRad::knots(self, tknots)
    }

    fn mults(&mut self, tmults: &mut [i32]) {
        BlendFuncEvolRad::mults(self, tmults)
    }

    fn section_d1(
        &mut self,
        p: &super::brep_blend_point::BlendPoint,
        poles: &mut [DVec3],
        d_poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        d_poles_2d: &mut [DVec2],
        weigths: &mut [f64],
        d_weigths: &mut [f64],
    ) -> bool {
        BlendFuncEvolRad::section_d1(
            self, p, poles, d_poles, poles_2d, d_poles_2d, weigths, d_weigths,
        )
    }

    fn section(
        &mut self,
        p: &super::brep_blend_point::BlendPoint,
        poles: &mut [DVec3],
        poles_2d: &mut [DVec2],
        weigths: &mut [f64],
    ) {
        BlendFuncEvolRad::section_simple(self, p, poles, poles_2d, weigths)
    }

    fn section_d2(
        &mut self,
        p: &super::brep_blend_point::BlendPoint,
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
        BlendFuncEvolRad::section_d2(
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
        BlendFuncEvolRad::resolution(self, ic_2d, tol, tol_u, tol_v)
    }
}

impl<'a> BlendFunction for BlendFuncEvolRad<'a> {
    fn point_on_s1(&self) -> DVec3 {
        BlendFuncEvolRad::point_on_s1(self)
    }

    fn point_on_s2(&self) -> DVec3 {
        BlendFuncEvolRad::point_on_s2(self)
    }

    fn is_tangency_point(&self) -> bool {
        BlendFuncEvolRad::is_tangency_point(self)
    }

    fn tangent_on_s1(&self) -> DVec3 {
        BlendFuncEvolRad::tangent_on_s1(self)
    }

    fn tangent_2d_on_s1(&self) -> DVec2 {
        BlendFuncEvolRad::tangent_2d_on_s1(self)
    }

    fn tangent_on_s2(&self) -> DVec3 {
        BlendFuncEvolRad::tangent_on_s2(self)
    }

    fn tangent_2d_on_s2(&self) -> DVec2 {
        BlendFuncEvolRad::tangent_2d_on_s2(self)
    }

    fn twist_on_s1(&self) -> bool {
        BlendFuncEvolRad::twist_on_s1(self)
    }

    fn twist_on_s2(&self) -> bool {
        BlendFuncEvolRad::twist_on_s2(self)
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
        BlendFuncEvolRad::tangent(self, u1, v1, u2, v2, tg_first, tg_last, norm_first, norm_last)
    }
}
