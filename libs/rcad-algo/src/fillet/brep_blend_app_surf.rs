//! OCCT BRepBlend AppSurf / AppFunc family (TKFillet/BRepBlend) —
//! AppBlend_Approx (TKGeomAlgo/AppBlend/AppBlend_Approx.hxx L32-93),
//! Approx_SweepFunction (TKGeomBase/Approx/Approx_SweepFunction.hxx L37-135),
//! BRepBlend_AppFuncRoot (BRepBlend_AppFuncRoot.hxx L37-165 +
//! BRepBlend_AppFuncRoot.cxx L26-433), BRepBlend_AppFunc (hxx L34-56 +
//! cxx L26-48), BRepBlend_AppFuncRst (cxx L26-48), BRepBlend_AppFuncRstRst
//! (cxx L26-48) and BRepBlend_AppSurface (hxx L37-125 + cxx L26-150 + .lxx).
//!
//! Architecture mappings: `occ::handle<T>` members -> owned values / shared
//! references; `math_Vector sol(1, N)` -> `Vec<f64>`; `NCollection_Array1`
//! out-parameters -> `&mut [T]` slices; the C++ inheritance chain
//! Approx_SweepFunction <- BRepBlend_AppFuncRoot <- AppFunc/AppFuncRst/
//! AppFuncRstRst maps to the [`ApproxSweepFunction`] trait implemented by the
//! three concrete structs over a shared [`BRepBlendAppFuncRoot`] core, with
//! the virtual Point/Vec overrides dispatched through `AppFuncKind`
//! (the OCCT hierarchy has exactly these three overrides).
//! `occ::handle<BRepBlend_Line>` (shared, mutated through InsertBefore) maps
//! to an owned clone inside the root: the ChFi3d callers already pass a
//! cloned line (chfi3d_builder_6b `the_walk.line().clone()`), so the
//! SearchPoint enrichment stays visible to the root's own binary searches,
//! which is the only consumer the algorithm relies on.
//!
//! PENDING (marked per plan 0.6, OCCT-named placeholders + failure path):
//! `Approx_SweepApproximation` (TKGeomAlgo/Approx) is not translated yet;
//! [`ApproxSweepApproximation`] stores the parameters and reports not-done,
//! so [`BRepBlendAppSurface::is_done`] returns false and the ChFi3d
//! CompleteData callers take their OCCT failure path.

use std::cell::RefCell;

use glam::{DVec2, DVec3};

use rcad_kernel::math::function_set_root::FunctionSetRoot;
use rcad_kernel::math::GeomAbsShape;

use super::brep_blend_function::BlendAppFunction;
use super::brep_blend_line::BRepBlendLine;
use super::brep_blend_point::BlendPoint;

/// OCCT AppBlend_Approx — deferred class for the bspline approximation of a
/// surface (AppBlend_Approx.hxx L32).  OCCT abstract class -> Rust trait.
pub trait AppBlendApprox {
    /// OCCT IsDone() (hxx L37).
    fn is_done(&self) -> bool;

    /// OCCT SurfShape(UDegree, VDegree, NbUPoles, NbVPoles, NbUKnots,
    /// NbVKnots) (hxx L39-44).
    fn surf_shape(
        &self,
        u_degree: &mut i32,
        v_degree: &mut i32,
        nb_u_poles: &mut i32,
        nb_v_poles: &mut i32,
        nb_u_knots: &mut i32,
        nb_v_knots: &mut i32,
    );

    /// OCCT Surface(TPoles, TWeights, TUKnots, TVKnots, TUMults, TVMults)
    /// (hxx L46-51).
    #[allow(clippy::too_many_arguments)]
    fn surface(
        &self,
        t_poles: &mut [Vec<DVec3>],
        t_weights: &mut [Vec<f64>],
        t_u_knots: &mut [f64],
        t_v_knots: &mut [f64],
        t_u_mults: &mut [i32],
        t_v_mults: &mut [i32],
    );

    /// OCCT UDegree() (hxx L53).
    fn u_degree(&self) -> i32;

    /// OCCT VDegree() (hxx L55).
    fn v_degree(&self) -> i32;

    /// OCCT SurfPoles() (hxx L57).
    fn surf_poles(&self) -> &[Vec<DVec3>];

    /// OCCT SurfWeights() (hxx L59).
    fn surf_weights(&self) -> &[Vec<f64>];

    /// OCCT SurfUKnots() (hxx L61).
    fn surf_u_knots(&self) -> &[f64];

    /// OCCT SurfVKnots() (hxx L63).
    fn surf_v_knots(&self) -> &[f64];

    /// OCCT SurfUMults() (hxx L65).
    fn surf_u_mults(&self) -> &[i32];

    /// OCCT SurfVMults() (hxx L67).
    fn surf_v_mults(&self) -> &[i32];

    /// OCCT NbCurves2d() (hxx L69).
    fn nb_curves2d(&self) -> i32;

    /// OCCT Curves2dShape(Degree, NbPoles, NbKnots) (hxx L71).
    fn curves2d_shape(&self, degree: &mut i32, nb_poles: &mut i32, nb_knots: &mut i32);

    /// OCCT Curve2d(Index, TPoles, TKnots, TMults) (hxx L73-76).
    fn curve2d(&self, index: i32, t_poles: &mut [DVec2], t_knots: &mut [f64], t_mults: &mut [i32]);

    /// OCCT Curves2dDegree() (hxx L78).
    fn curves2d_degree(&self) -> i32;

    /// OCCT Curve2dPoles(Index) (hxx L80-81).
    fn curve2d_poles(&self, index: i32) -> &[DVec2];

    /// OCCT Curves2dKnots() (hxx L83).
    fn curves2d_knots(&self) -> &[f64];

    /// OCCT Curves2dMults() (hxx L85).
    fn curves2d_mults(&self) -> &[i32];

    /// OCCT TolReached(Tol3d, Tol2d) (hxx L87).
    fn tol_reached(&self, tol3d: &mut f64, tol2d: &mut f64);

    /// OCCT TolCurveOnSurf(Index) (hxx L89).
    fn tol_curve_on_surf(&self, index: i32) -> f64;
}

/// OCCT Approx_SweepFunction — the function used by SweepApproximation to
/// perform sweeping application (Approx_SweepFunction.hxx L37).  OCCT
/// abstract class -> Rust trait (the D1/D2/Resolution/BarycentreOfSurf/
/// MaximalSection/GetMinimalWeight OCCT default implementations live in
/// AppFuncRoot and are forwarded by it).
///
/// Architecture note: the OCCT entry points are declared non-const, but the
/// consumers pass the function as `const handle<Approx_SweepFunction>&`
/// (the pointed object stays mutable through the handle); the Rust trait
/// therefore takes `&self` and the concrete adapters keep their mutable
/// state behind a RefCell (see BRepBlendAppFuncRoot).
pub trait ApproxSweepFunction {
    /// OCCT D0(Param, First, Last, Poles, Poles2d, Weigths) (hxx L43-49).
    #[allow(clippy::too_many_arguments)]
    fn d0(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        weigths: &mut [f64],
    ) -> bool;

    /// OCCT D1(Param, First, Last, Poles, DPoles, Poles2d, DPoles2d,
    /// Weigths, DWeigths) (hxx L55-63).
    #[allow(clippy::too_many_arguments)]
    fn d1(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool;

    /// OCCT D2(Param, First, Last, Poles, DPoles, D2Poles, Poles2d,
    /// DPoles2d, D2Poles2d, Weigths, DWeigths, D2Weigths) (hxx L67-78).
    #[allow(clippy::too_many_arguments)]
    fn d2(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        d2poles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
        d2weigths: &mut [f64],
    ) -> bool;

    /// OCCT Nb2dCurves() (hxx L81).
    fn nb2dcurves(&self) -> i32;

    /// OCCT SectionShape(NbPoles, NbKnots, Degree) (hxx L84).
    fn section_shape(&self, nb_poles: &mut i32, nb_knots: &mut i32, degree: &mut i32);

    /// OCCT Knots(TKnots) (hxx L87).
    fn knots(&self, tknots: &mut [f64]);

    /// OCCT Mults(TMults) (hxx L90).
    fn mults(&self, tmults: &mut [i32]);

    /// OCCT IsRational() (hxx L93).
    fn is_rational(&self) -> bool;

    /// OCCT NbIntervals(S) (hxx L97).
    fn nb_intervals(&self, s: GeomAbsShape) -> i32;

    /// OCCT Intervals(T, S) (hxx L104-105).
    fn intervals(&self, t: &mut [f64], s: GeomAbsShape);

    /// OCCT SetInterval(First, Last) (hxx L111).
    fn set_interval(&self, first: f64, last: f64);

    /// OCCT Resolution(Index, Tol, TolU, TolV) (hxx L116-119).
    fn resolution(&self, index: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64);

    /// OCCT GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d) (hxx L126-129).
    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64, tol3d: &mut [f64]);

    /// OCCT SetTolerance(Tol3d, Tol2d) (hxx L133).
    fn set_tolerance(&self, tol3d: f64, tol2d: f64);

    /// OCCT BarycentreOfSurf() (hxx L138).
    fn barycentre_of_surf(&self) -> DVec3;

    /// OCCT MaximalSection() (hxx L145).
    fn maximal_section(&self) -> f64;

    /// OCCT GetMinimalWeight(Weigths) (hxx L150).
    fn get_minimal_weight(&self, weigths: &mut [f64]);
}

/// Which OCCT Point/Vec override the concrete adapter carries (the virtual
/// table of the AppFuncRoot hierarchy; OCCT has exactly these three).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AppFuncKind {
    /// OCCT BRepBlend_AppFunc (surface / surface).
    Func,
    /// OCCT BRepBlend_AppFuncRst (surface / restriction).
    FuncRst,
    /// OCCT BRepBlend_AppFuncRstRst (restriction / restriction).
    FuncRstRst,
}

/// Mutable core state of BRepBlend_AppFuncRoot (BRepBlend_AppFuncRoot.hxx
/// L155-164 private fields).
struct RootState<'a> {
    /// OCCT: occ::handle<BRepBlend_Line> myLine (owned clone; see the module
    /// architecture note).
    my_line: BRepBlendLine,
    /// OCCT: void* myFunc (a Blend_AppFunction&).
    my_func: &'a mut dyn BlendAppFunction,
    /// OCCT: math_Vector myTolerance.
    my_tolerance: Vec<f64>,
    /// OCCT: Blend_Point myPnt.
    my_pnt: BlendPoint,
    /// OCCT: gp_Pnt myBary.
    my_bary: DVec3,
    /// OCCT: math_Vector X1.
    x1: Vec<f64>,
    /// OCCT: math_Vector X2.
    x2: Vec<f64>,
    /// OCCT: math_Vector XInit.
    x_init: Vec<f64>,
    /// OCCT: math_Vector Sol.
    sol: Vec<f64>,
    /// OCCT: std::unique_ptr<math_FunctionSetRoot> mySolver.
    my_solver: FunctionSetRoot,
}

/// OCCT BRepBlend_AppFuncRoot — root class for the function-to-approximation
/// adapters (BRepBlend_AppFuncRoot.hxx L37).  The mutable state lives behind
/// a RefCell because the OCCT Approx_SweepFunction entry points are
/// non-const while the ChFi3d callers hold the adapter through a shared
/// reference when handing it to BRepBlend_AppSurface.
pub struct BRepBlendAppFuncRoot<'a> {
    state: RefCell<RootState<'a>>,
    /// The Point/Vec override this instance carries.
    kind: AppFuncKind,
}

impl<'a> BRepBlendAppFuncRoot<'a> {
    /// OCCT BRepBlend_AppFuncRoot(Line, Func, Tol3d, Tol2d)
    /// (BRepBlend_AppFuncRoot.cxx L26-80).
    pub fn new(
        line: &'a BRepBlendLine,
        func: &'a mut dyn BlendAppFunction,
        tol3d: f64,
        tol2d: f64,
    ) -> Self {
        // OCCT: myTolerance(1, Func.NbVariables()), X1/X2/XInit/Sol likewise.
        let dim = func.nb_variables();
        let mut my_tolerance = vec![0.0; dim];
        let x1 = vec![0.0; dim];
        let x2 = vec![0.0; dim];
        let x_init = vec![0.0; dim];
        let sol = vec![0.0; dim];

        //  Tolerances
        // OCCT: Func.GetTolerance(myTolerance, Tol3d);
        func.get_tolerance(&mut my_tolerance, tol3d);
        // OCCT: for (ii = 1; ii <= dim; ii++) if (myTolerance(ii) > Tol2d) ...
        for t in my_tolerance.iter_mut() {
            if *t > tol2d {
                *t = tol2d;
            }
        }

        //  Tables
        // OCCT: Func.GetShape(NbPoles, NbKnots, Degree, NbPoles2d);
        let mut nb_poles = 0i32;
        let mut nb_knots = 0i32;
        let mut degree = 0i32;
        let mut nb_poles_2d = 0i32;
        func.get_shape(&mut nb_poles, &mut nb_knots, &mut degree, &mut nb_poles_2d);

        // Calculation of BaryCentre (rational case).
        let mut my_bary = DVec3::ZERO;
        if func.is_rational() {
            let mut xmax: f64 = -1.0e100;
            let mut xmin: f64 = 1.0e100;
            let mut ymax: f64 = -1.0e100;
            let mut ymin: f64 = 1.0e100;
            let mut zmax: f64 = -1.0e100;
            let mut zmin: f64 = 1.0e100;
            let mut p;
            for ii in 1..=line.nb_points() {
                p = line.point(ii).clone();
                xmax = xmax.max(p.point_on_s1().x.max(p.point_on_s2().x));
                xmin = xmin.min(p.point_on_s1().x.min(p.point_on_s2().x));
                ymax = ymax.max(p.point_on_s1().y.max(p.point_on_s2().y));
                ymin = ymin.min(p.point_on_s1().y.min(p.point_on_s2().y));
                zmax = zmax.max(p.point_on_s1().z.max(p.point_on_s2().z));
                zmin = zmin.min(p.point_on_s1().z.min(p.point_on_s2().z));

                my_bary = DVec3::new(
                    (xmax + xmin) / 2.0,
                    (ymax + ymin) / 2.0,
                    (zmax + zmin) / 2.0,
                );
            }
        } else {
            my_bary = DVec3::ZERO;
        }

        // OCCT: mySolver = std::make_unique<math_FunctionSetRoot>(Func,
        // myTolerance, 30);
        let fswd: &dyn rcad_kernel::math::function_set_root::FunctionSetWithDerivatives = &*func;
        let my_solver = FunctionSetRoot::new(fswd, &my_tolerance, 30);

        BRepBlendAppFuncRoot {
            state: RefCell::new(RootState {
                my_line: line.clone(),
                my_func: func,
                my_tolerance,
                my_pnt: BlendPoint::new(),
                my_bary,
                x1,
                x2,
                x_init,
                sol,
                my_solver,
            }),
            // The concrete adapter (AppFunc/AppFuncRst/AppFuncRstRst) sets
            // the override after construction; default to the S/S table.
            kind: AppFuncKind::Func,
        }
    }

    /// Bind the Point/Vec override (the C++ vtable binding done by the
    /// concrete derived constructors).
    fn set_kind(&mut self, kind: AppFuncKind) {
        self.kind = kind;
    }

    /// OCCT D0(Param, First, Last, Poles, Poles2d, Weigths)
    /// (BRepBlend_AppFuncRoot.cxx L91-107).
    #[allow(clippy::too_many_arguments)]
    pub fn d0(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
        poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        weigths: &mut [f64],
    ) -> bool {
        let mut st = self.state.borrow_mut();
        // OCCT: Ok = SearchPoint(*Func, Param, myPnt);
        let mut my_pnt = st.my_pnt.clone();
        let ok = Self::search_point(self.kind, &mut st, param, &mut my_pnt);
        st.my_pnt = my_pnt;
        if ok {
            // OCCT: (*Func).Section(myPnt, Poles, Poles2d, Weigths);
            let pnt = st.my_pnt.clone();
            let func = &mut *st.my_func;
            func.section(&pnt, poles, poles2d, weigths);
        }
        ok
    }

    /// OCCT D1(Param, First, Last, Poles, DPoles, Poles2d, DPoles2d,
    /// Weigths, DWeigths) (BRepBlend_AppFuncRoot.cxx L114-135).
    #[allow(clippy::too_many_arguments)]
    pub fn d1(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool {
        let mut st = self.state.borrow_mut();
        let mut my_pnt = st.my_pnt.clone();
        let mut ok = Self::search_point(self.kind, &mut st, param, &mut my_pnt);
        st.my_pnt = my_pnt;
        if ok {
            // OCCT: Ok = (*Func).Section(myPnt, Poles, DPoles, Poles2d,
            // DPoles2d, Weigths, DWeigths);
            let pnt = st.my_pnt.clone();
            let func = &mut *st.my_func;
            ok = func.section_d1(&pnt, poles, dpoles, poles2d, dpoles2d, weigths, dweigths);
        }
        ok
    }

    /// OCCT D2(Param, First, Last, Poles, DPoles, D2Poles, Poles2d,
    /// DPoles2d, D2Poles2d, Weigths, DWeigths, D2Weigths)
    /// (BRepBlend_AppFuncRoot.cxx L143-174).
    #[allow(clippy::too_many_arguments)]
    pub fn d2(
        &self,
        param: f64,
        _first: f64,
        _last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        d2poles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
        d2weigths: &mut [f64],
    ) -> bool {
        let mut st = self.state.borrow_mut();
        let mut my_pnt = st.my_pnt.clone();
        let mut ok = Self::search_point(self.kind, &mut st, param, &mut my_pnt);
        st.my_pnt = my_pnt;
        if ok {
            // OCCT: Ok = (*Func).Section(myPnt, Poles, DPoles, D2Poles,
            // Poles2d, DPoles2d, D2Poles2d, Weigths, DWeigths, D2Weigths);
            let pnt = st.my_pnt.clone();
            let func = &mut *st.my_func;
            ok = func.section_d2(
                &pnt,
                poles,
                dpoles,
                d2poles,
                poles2d,
                dpoles2d,
                d2poles2d,
                weigths,
                dweigths,
                d2weigths,
            );
        }
        ok
    }

    /// OCCT Nb2dCurves() (BRepBlend_AppFuncRoot.cxx L176-182).
    pub fn nb_2d_curves(&self) -> i32 {
        let mut st = self.state.borrow_mut();
        let func = &mut *st.my_func;
        let mut i = 0i32;
        let mut j = 0i32;
        let mut k = 0i32;
        let mut nbpol2d = 0i32;
        func.get_shape(&mut i, &mut j, &mut k, &mut nbpol2d);
        nbpol2d
    }

    /// OCCT SectionShape(NbPoles, NbKnots, Degree)
    /// (BRepBlend_AppFuncRoot.cxx L184-189).
    pub fn section_shape(&self, nb_poles: &mut i32, nb_knots: &mut i32, degree: &mut i32) {
        let mut st = self.state.borrow_mut();
        let func = &mut *st.my_func;
        let mut ii = 0i32;
        func.get_shape(nb_poles, nb_knots, degree, &mut ii);
    }

    /// OCCT Knots(TKnots) (BRepBlend_AppFuncRoot.cxx L191-195).
    pub fn knots(&self, tknots: &mut [f64]) {
        let mut st = self.state.borrow_mut();
        st.my_func.knots(tknots);
    }

    /// OCCT Mults(TMults) (BRepBlend_AppFuncRoot.cxx L197-201).
    pub fn mults(&self, tmults: &mut [i32]) {
        let mut st = self.state.borrow_mut();
        st.my_func.mults(tmults);
    }

    /// OCCT IsRational() (BRepBlend_AppFuncRoot.cxx L203-207).
    pub fn is_rational(&self) -> bool {
        let st = self.state.borrow();
        st.my_func.is_rational()
    }

    /// OCCT NbIntervals(S) (BRepBlend_AppFuncRoot.cxx L209-213).
    pub fn nb_intervals(&self, s: GeomAbsShape) -> i32 {
        let st = self.state.borrow();
        st.my_func.nb_intervals(s) as i32
    }

    /// OCCT Intervals(T, S) (BRepBlend_AppFuncRoot.cxx L215-219).
    pub fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        let st = self.state.borrow();
        st.my_func.intervals(t, s);
    }

    /// OCCT SetInterval(First, Last) (BRepBlend_AppFuncRoot.cxx L221-225)
    /// — Func->Set(First, Last).
    pub fn set_interval(&self, first: f64, last: f64) {
        let mut st = self.state.borrow_mut();
        st.my_func.set_interval(first, last);
    }

    /// OCCT Resolution(Index, Tol, TolU, TolV)
    /// (BRepBlend_AppFuncRoot.cxx L227-234).
    pub fn resolution(&self, index: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        let st = self.state.borrow();
        st.my_func.resolution(index, tol, tol_u, tol_v);
    }

    /// OCCT GetTolerance(BoundTol, SurfTol, AngleTol, Tol3d)
    /// (BRepBlend_AppFuncRoot.cxx L236-250).
    pub fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64, tol3d: &mut [f64]) {
        let st = self.state.borrow();
        let len = tol3d.len();
        let mut v3d = vec![0.0; len];
        let mut v1d = vec![0.0; len];
        st.my_func
            .get_approx_tolerance(bound_tol, surf_tol, angle_tol, &mut v3d, &mut v1d);
        // OCCT: for (ii = 1; ii <= Tol3d.Length(); ii++) Tol3d(ii) = V3d(ii);
        tol3d.copy_from_slice(&v3d);
    }

    /// OCCT SetTolerance(Tol3d, Tol2d) (BRepBlend_AppFuncRoot.cxx L252-265).
    pub fn set_tolerance(&self, tol3d: f64, tol2d: f64) {
        let mut st = self.state.borrow_mut();
        let dim = st.my_func.nb_variables();
        let mut my_tolerance = std::mem::take(&mut st.my_tolerance);
        st.my_func.get_tolerance(&mut my_tolerance, tol3d);
        for ii in 0..dim {
            if my_tolerance[ii] > tol2d {
                my_tolerance[ii] = tol2d;
            }
        }
        st.my_tolerance = my_tolerance;
        // OCCT: mySolver->SetTolerance(myTolerance);
        let tolerance = st.my_tolerance.clone();
        st.my_solver.set_tolerance(&tolerance);
    }

    /// OCCT BarycentreOfSurf() (BRepBlend_AppFuncRoot.cxx L267-270).
    pub fn barycentre_of_surf(&self) -> DVec3 {
        let st = self.state.borrow();
        st.my_bary
    }

    /// OCCT MaximalSection() (BRepBlend_AppFuncRoot.cxx L272-276) —
    /// Func->GetSectionSize().
    pub fn maximal_section(&self) -> f64 {
        let st = self.state.borrow();
        st.my_func.get_section_size()
    }

    /// OCCT GetMinimalWeight(Weigths) (BRepBlend_AppFuncRoot.cxx L278-282).
    pub fn get_minimal_weight(&self, weigths: &mut [f64]) {
        let st = self.state.borrow();
        st.my_func.get_minimal_weight(weigths);
    }

    /// OCCT Point(Func, Param, Sol, Pnt) — virtual dispatch over the three
    /// derived overrides.
    fn point_dispatch(
        kind: AppFuncKind,
        func: &dyn BlendAppFunction,
        param: f64,
        the_sol: &[f64],
        pnt: &mut BlendPoint,
    ) {
        match kind {
            // OCCT BRepBlend_AppFunc::Point (BRepBlend_AppFunc.cxx L36-42):
            // Pnt.SetValue(Func.Pnt1(), Func.Pnt2(), Param, theSol(1..4)).
            AppFuncKind::Func => {
                pnt.set_value_on_2_surfaces(
                    func.pnt1(),
                    func.pnt2(),
                    param,
                    the_sol[0],
                    the_sol[1],
                    the_sol[2],
                    the_sol[3],
                );
            }
            // OCCT BRepBlend_AppFuncRst::Point (BRepBlend_AppFuncRst.cxx
            // L36-42):
            // Pnt.SetValue(Func.Pnt1(), Func.Pnt2(), Param, theSol(1..3)).
            AppFuncKind::FuncRst => {
                pnt.set_value_on_surface_curve(
                    func.pnt1(),
                    func.pnt2(),
                    param,
                    the_sol[0],
                    the_sol[1],
                    the_sol[2],
                );
            }
            // OCCT BRepBlend_AppFuncRstRst::Point (BRepBlend_AppFuncRstRst.cxx
            // L36-41):
            // Pnt.SetValue(Func.Pnt1(), Func.Pnt2(), Param, theSol(1..2)).
            AppFuncKind::FuncRstRst => {
                pnt.set_value_on_2_curves(func.pnt1(), func.pnt2(), param, the_sol[0], the_sol[1]);
            }
        }
    }

    /// OCCT Vec(Sol, Pnt) — virtual dispatch over the three derived
    /// overrides (fills the solution vector from the point parameters).
    fn vec_dispatch(kind: AppFuncKind, the_sol: &mut [f64], pnt: &BlendPoint) {
        match kind {
            // OCCT BRepBlend_AppFunc::Vec (BRepBlend_AppFunc.cxx L44-48):
            // Pnt.ParametersOnS1(theSol(1), theSol(2));
            // Pnt.ParametersOnS2(theSol(3), theSol(4)).
            AppFuncKind::Func => {
                let (u1, v1) = pnt.parameters_on_s1();
                the_sol[0] = u1;
                the_sol[1] = v1;
                let (u2, v2) = pnt.parameters_on_s2();
                the_sol[2] = u2;
                the_sol[3] = v2;
            }
            // OCCT BRepBlend_AppFuncRst::Vec (BRepBlend_AppFuncRst.cxx
            // L44-48):
            // Pnt.ParametersOnS(theSol(1), theSol(2));
            // theSol(3) = Pnt.ParameterOnC().
            AppFuncKind::FuncRst => {
                let (u, v) = pnt.parameters_on_s();
                the_sol[0] = u;
                the_sol[1] = v;
                the_sol[2] = pnt.parameter_on_c();
            }
            // OCCT BRepBlend_AppFuncRstRst::Vec (BRepBlend_AppFuncRstRst.cxx
            // L43-47):
            // theSol(1) = Pnt.ParameterOnC1();
            // theSol(2) = Pnt.ParameterOnC2().
            AppFuncKind::FuncRstRst => {
                the_sol[0] = pnt.parameter_on_c1();
                the_sol[1] = pnt.parameter_on_c2();
            }
        }
    }

    /// OCCT SearchPoint(Func, Param, Pnt) (BRepBlend_AppFuncRoot.cxx
    /// L300-375) — find the solution point with parameter Param (on the 2
    /// supports).
    fn search_point(
        kind: AppFuncKind,
        st: &mut RootState<'a>,
        param: f64,
        pnt: &mut BlendPoint,
    ) -> bool {
        let dim = st.my_func.nb_variables();
        // (1) Find a point of init
        let index;
        let i2 = st.my_line.nb_points();

        //  (1.a) It is checked if it is inside
        if param < st.my_line.point(1).parameter() {
            return false;
        }
        if param > st.my_line.point(i2).parameter() {
            return false;
        }

        //  (1.b) Find the interval
        let mut param_index = 0i32;
        let trouve = Self::search_location(st, param, 1, i2, &mut param_index);
        index = param_index;

        //  (1.c) If the point is already calculated it is returned
        if trouve {
            *pnt = st.my_line.point(index).clone();
            // OCCT: Vec(XInit, Pnt);
            let mut x_init = std::mem::take(&mut st.x_init);
            Self::vec_dispatch(kind, &mut x_init, pnt);
            st.x_init = x_init;
        } else {
            //  (1.d) Initialisation by linear interpolation
            *pnt = st.my_line.point(index).clone();
            // OCCT: Vec(X1, Pnt); t1 = Pnt.Parameter();
            let mut x1 = std::mem::take(&mut st.x1);
            Self::vec_dispatch(kind, &mut x1, pnt);
            st.x1 = x1;
            let t1 = pnt.parameter();

            *pnt = st.my_line.point(index + 1).clone();
            // OCCT: Vec(X2, Pnt); t2 = Pnt.Parameter();
            let mut x2 = std::mem::take(&mut st.x2);
            Self::vec_dispatch(kind, &mut x2, pnt);
            st.x2 = x2;
            let t2 = pnt.parameter();

            let parammt1 = (param - t1) / (t2 - t1);
            let t2mparam = (t2 - param) / (t2 - t1);
            for i in 0..dim {
                st.x_init[i] = st.x2[i] * parammt1 + st.x1[i] * t2mparam;
            }
        }

        // (2) Calculation of the solution ------------------------
        let func = &mut *st.my_func;
        func.set_param(param);
        // OCCT: Func.GetBounds(X1, X2);
        func.get_bounds(&mut st.x1, &mut st.x2);
        // OCCT: mySolver->Perform(Func, XInit, X1, X2);
        st.my_solver.perform(func, &st.x_init, &st.x1, &st.x2, false);

        if !st.my_solver.is_done() {
            return false;
        }
        // OCCT: mySolver->Root(Sol);
        let root = st.my_solver.root();
        st.sol.copy_from_slice(&root);

        // (3) Storage of the point
        {
            let func = &mut *st.my_func;
            let sol = st.sol.clone();
            Self::point_dispatch(kind, func, param, &sol, pnt);
        }

        // (4) Insertion of the point if the calculation seems long.
        if !trouve && st.my_solver.nb_iterations() > 3 {
            // OCCT: myLine->InsertBefore(Index + 1, Pnt);
            st.my_line.insert_before(index + 1, pnt.clone());
        }
        true
    }

    /// OCCT SearchLocation(Param, FirstIndex, LastIndex, ParamIndex)
    /// (BRepBlend_AppFuncRoot.cxx L387-433) — binary search of the line for
    /// the parametric interval containing Param.
    fn search_location(
        st: &RootState<'a>,
        param: f64,
        first_index: i32,
        last_index: i32,
        param_index: &mut i32,
    ) -> bool {
        let mut ideb = first_index;
        let mut ifin = last_index;
        let mut idemi;
        let mut valeur;

        valeur = st.my_line.point(ideb).parameter();
        if param == valeur {
            *param_index = ideb;
            return true;
        }

        valeur = st.my_line.point(ifin).parameter();
        if param == valeur {
            *param_index = ifin;
            return true;
        }

        while ideb + 1 != ifin {
            idemi = (ideb + ifin) / 2;
            valeur = st.my_line.point(idemi).parameter();
            if valeur < param {
                ideb = idemi;
            } else if valeur > param {
                ifin = idemi;
            } else {
                *param_index = idemi;
                return true;
            }
        }

        *param_index = ideb;
        false
    }
}

/// OCCT BRepBlend_AppFunc — adapter binding a Blend_Function
/// (surface / surface blending) to the approximation driver
/// (BRepBlend_AppFunc.hxx L34; BRepBlend_AppFunc.cxx L26-33).
pub struct BRepBlendAppFunc<'a> {
    pub(crate) root: BRepBlendAppFuncRoot<'a>,
}

impl<'a> BRepBlendAppFunc<'a> {
    /// OCCT BRepBlend_AppFunc(Line, Func, Tol3d, Tol2d)
    /// (BRepBlend_AppFunc.cxx L26-33).
    pub fn new(
        line: &'a BRepBlendLine,
        func: &'a mut dyn BlendAppFunction,
        tol3d: f64,
        tol2d: f64,
    ) -> Self {
        let mut root = BRepBlendAppFuncRoot::new(line, func, tol3d, tol2d);
        root.set_kind(AppFuncKind::Func);
        BRepBlendAppFunc { root }
    }
}

impl ApproxSweepFunction for BRepBlendAppFunc<'_> {
    fn d0(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        weigths: &mut [f64],
    ) -> bool {
        self.root.d0(param, first, last, poles, poles2d, weigths)
    }

    fn d1(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool {
        self.root
            .d1(param, first, last, poles, dpoles, poles2d, dpoles2d, weigths, dweigths)
    }

    fn d2(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        d2poles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
        d2weigths: &mut [f64],
    ) -> bool {
        self.root.d2(
            param,
            first,
            last,
            poles,
            dpoles,
            d2poles,
            poles2d,
            dpoles2d,
            d2poles2d,
            weigths,
            dweigths,
            d2weigths,
        )
    }

    fn nb2dcurves(&self) -> i32 {
        self.root.nb_2d_curves()
    }

    fn section_shape(&self, nb_poles: &mut i32, nb_knots: &mut i32, degree: &mut i32) {
        self.root.section_shape(nb_poles, nb_knots, degree)
    }

    fn knots(&self, tknots: &mut [f64]) {
        self.root.knots(tknots)
    }

    fn mults(&self, tmults: &mut [i32]) {
        self.root.mults(tmults)
    }

    fn is_rational(&self) -> bool {
        self.root.is_rational()
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> i32 {
        self.root.nb_intervals(s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        self.root.intervals(t, s)
    }

    fn set_interval(&self, first: f64, last: f64) {
        self.root.set_interval(first, last)
    }

    fn resolution(&self, index: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.root.resolution(index, tol, tol_u, tol_v)
    }

    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64, tol3d: &mut [f64]) {
        self.root.get_tolerance(bound_tol, surf_tol, angle_tol, tol3d)
    }

    fn set_tolerance(&self, tol3d: f64, tol2d: f64) {
        self.root.set_tolerance(tol3d, tol2d)
    }

    fn barycentre_of_surf(&self) -> DVec3 {
        self.root.barycentre_of_surf()
    }

    fn maximal_section(&self) -> f64 {
        self.root.maximal_section()
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        self.root.get_minimal_weight(weigths)
    }
}

/// OCCT BRepBlend_AppFuncRst — adapter binding a Blend_SurfRstFunction
/// (surface / restriction blending) to the approximation driver
/// (BRepBlend_AppFuncRst.hxx L34; BRepBlend_AppFuncRst.cxx L26-34).
pub struct BRepBlendAppFuncRst<'a> {
    pub(crate) root: BRepBlendAppFuncRoot<'a>,
}

impl<'a> BRepBlendAppFuncRst<'a> {
    /// OCCT BRepBlend_AppFuncRst(Line, Func, Tol3d, Tol2d)
    /// (BRepBlend_AppFuncRst.cxx L26-34).
    pub fn new(
        line: &'a BRepBlendLine,
        func: &'a mut dyn BlendAppFunction,
        tol3d: f64,
        tol2d: f64,
    ) -> Self {
        let mut root = BRepBlendAppFuncRoot::new(line, func, tol3d, tol2d);
        root.set_kind(AppFuncKind::FuncRst);
        BRepBlendAppFuncRst { root }
    }
}

impl ApproxSweepFunction for BRepBlendAppFuncRst<'_> {
    fn d0(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        weigths: &mut [f64],
    ) -> bool {
        self.root.d0(param, first, last, poles, poles2d, weigths)
    }

    fn d1(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool {
        self.root
            .d1(param, first, last, poles, dpoles, poles2d, dpoles2d, weigths, dweigths)
    }

    fn d2(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        d2poles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
        d2weigths: &mut [f64],
    ) -> bool {
        self.root.d2(
            param,
            first,
            last,
            poles,
            dpoles,
            d2poles,
            poles2d,
            dpoles2d,
            d2poles2d,
            weigths,
            dweigths,
            d2weigths,
        )
    }

    fn nb2dcurves(&self) -> i32 {
        self.root.nb_2d_curves()
    }

    fn section_shape(&self, nb_poles: &mut i32, nb_knots: &mut i32, degree: &mut i32) {
        self.root.section_shape(nb_poles, nb_knots, degree)
    }

    fn knots(&self, tknots: &mut [f64]) {
        self.root.knots(tknots)
    }

    fn mults(&self, tmults: &mut [i32]) {
        self.root.mults(tmults)
    }

    fn is_rational(&self) -> bool {
        self.root.is_rational()
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> i32 {
        self.root.nb_intervals(s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        self.root.intervals(t, s)
    }

    fn set_interval(&self, first: f64, last: f64) {
        self.root.set_interval(first, last)
    }

    fn resolution(&self, index: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.root.resolution(index, tol, tol_u, tol_v)
    }

    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64, tol3d: &mut [f64]) {
        self.root.get_tolerance(bound_tol, surf_tol, angle_tol, tol3d)
    }

    fn set_tolerance(&self, tol3d: f64, tol2d: f64) {
        self.root.set_tolerance(tol3d, tol2d)
    }

    fn barycentre_of_surf(&self) -> DVec3 {
        self.root.barycentre_of_surf()
    }

    fn maximal_section(&self) -> f64 {
        self.root.maximal_section()
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        self.root.get_minimal_weight(weigths)
    }
}

/// OCCT BRepBlend_AppFuncRstRst — adapter binding a Blend_RstRstFunction
/// (restriction / restriction blending) to the approximation driver
/// (BRepBlend_AppFuncRstRst.hxx L34; BRepBlend_AppFuncRstRst.cxx L26-34).
pub struct BRepBlendAppFuncRstRst<'a> {
    pub(crate) root: BRepBlendAppFuncRoot<'a>,
}

impl<'a> BRepBlendAppFuncRstRst<'a> {
    /// OCCT BRepBlend_AppFuncRstRst(Line, Func, Tol3d, Tol2d)
    /// (BRepBlend_AppFuncRstRst.cxx L26-34).
    pub fn new(
        line: &'a BRepBlendLine,
        func: &'a mut dyn BlendAppFunction,
        tol3d: f64,
        tol2d: f64,
    ) -> Self {
        let mut root = BRepBlendAppFuncRoot::new(line, func, tol3d, tol2d);
        root.set_kind(AppFuncKind::FuncRstRst);
        BRepBlendAppFuncRstRst { root }
    }
}

impl ApproxSweepFunction for BRepBlendAppFuncRstRst<'_> {
    fn d0(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        weigths: &mut [f64],
    ) -> bool {
        self.root.d0(param, first, last, poles, poles2d, weigths)
    }

    fn d1(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool {
        self.root
            .d1(param, first, last, poles, dpoles, poles2d, dpoles2d, weigths, dweigths)
    }

    fn d2(
        &self,
        param: f64,
        first: f64,
        last: f64,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        d2poles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        d2poles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
        d2weigths: &mut [f64],
    ) -> bool {
        self.root.d2(
            param,
            first,
            last,
            poles,
            dpoles,
            d2poles,
            poles2d,
            dpoles2d,
            d2poles2d,
            weigths,
            dweigths,
            d2weigths,
        )
    }

    fn nb2dcurves(&self) -> i32 {
        self.root.nb_2d_curves()
    }

    fn section_shape(&self, nb_poles: &mut i32, nb_knots: &mut i32, degree: &mut i32) {
        self.root.section_shape(nb_poles, nb_knots, degree)
    }

    fn knots(&self, tknots: &mut [f64]) {
        self.root.knots(tknots)
    }

    fn mults(&self, tmults: &mut [i32]) {
        self.root.mults(tmults)
    }

    fn is_rational(&self) -> bool {
        self.root.is_rational()
    }

    fn nb_intervals(&self, s: GeomAbsShape) -> i32 {
        self.root.nb_intervals(s)
    }

    fn intervals(&self, t: &mut [f64], s: GeomAbsShape) {
        self.root.intervals(t, s)
    }

    fn set_interval(&self, first: f64, last: f64) {
        self.root.set_interval(first, last)
    }

    fn resolution(&self, index: i32, tol: f64, tol_u: &mut f64, tol_v: &mut f64) {
        self.root.resolution(index, tol, tol_u, tol_v)
    }

    fn get_tolerance(&self, bound_tol: f64, surf_tol: f64, angle_tol: f64, tol3d: &mut [f64]) {
        self.root.get_tolerance(bound_tol, surf_tol, angle_tol, tol3d)
    }

    fn set_tolerance(&self, tol3d: f64, tol2d: f64) {
        self.root.set_tolerance(tol3d, tol2d)
    }

    fn barycentre_of_surf(&self) -> DVec3 {
        self.root.barycentre_of_surf()
    }

    fn maximal_section(&self) -> f64 {
        self.root.maximal_section()
    }

    fn get_minimal_weight(&self, weigths: &mut [f64]) {
        self.root.get_minimal_weight(weigths)
    }
}

/// OCCT Approx_SweepApproximation (TKGeomAlgo/Approx/
/// Approx_SweepApproximation.cxx) — PENDING TRANSLATION (plan 0.6).  The
/// placeholder stores the function reference and the Perform parameters and
/// reports not-done, which drives the OCCT failure path of the consumers
/// (BRepBlend_AppSurface::IsDone -> ChFi3d CompleteData returns false).
pub struct ApproxSweepApproximation<'a> {
    funct: &'a dyn ApproxSweepFunction,
    first: f64,
    last: f64,
    tol3d: f64,
    tol2d: f64,
    tol_angular: f64,
    continuity: crate::geomalgo::gtests_stubs::GeomAbsShape,
    degmax: i32,
    segmax: i32,
    done: bool,
}

impl<'a> ApproxSweepApproximation<'a> {
    /// OCCT Approx_SweepApproximation(Funct) — the construction part; the
    /// Perform(First, Last, Tol3d, Tol3d, Tol2d, TolAngular, Continuity,
    /// Degmax, Segmax) sweep driver is the pending boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn new(funct: &'a dyn ApproxSweepFunction) -> Self {
        ApproxSweepApproximation {
            funct,
            first: 0.0,
            last: 0.0,
            tol3d: 0.0,
            tol2d: 0.0,
            tol_angular: 0.0,
            continuity: crate::geomalgo::gtests_stubs::GeomAbsShape::C0,
            degmax: 11,
            segmax: 50,
            done: false,
        }
    }

    /// OCCT Approx_SweepApproximation::Perform — PENDING: the sweep
    /// approximation driver (Approx_SweepApproximation.cxx L100+).  The
    /// parameters are recorded; the algorithm reports not-done.
    #[allow(clippy::too_many_arguments)]
    pub fn perform(
        &mut self,
        first: f64,
        last: f64,
        tol3d: f64,
        boundary_tol: f64,
        tol2d: f64,
        tol_angular: f64,
        continuity: crate::geomalgo::gtests_stubs::GeomAbsShape,
        degmax: i32,
        segmax: i32,
    ) {
        self.first = first;
        self.last = last;
        self.tol3d = tol3d;
        self.tol2d = tol2d;
        self.tol_angular = tol_angular;
        self.continuity = continuity;
        self.degmax = degmax;
        self.segmax = segmax;
        let _ = boundary_tol;
        self.done = false;
    }

    /// OCCT IsDone().
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// The approximated function (read-only access for the delegated
    /// accessors of BRepBlend_AppSurface).
    pub(crate) fn funct(&self) -> &dyn ApproxSweepFunction {
        self.funct
    }
}

/// OCCT BRepBlend_AppSurface — used to approximate the blending surfaces
/// (BRepBlend_AppSurface.hxx L37; BRepBlend_AppSurface.cxx L26-150; the
/// inline accessors of BRepBlend_AppSurface.lxx delegate to approx).
pub struct BRepBlendAppSurface<'a> {
    approx: ApproxSweepApproximation<'a>,
}

impl<'a> BRepBlendAppSurface<'a> {
    /// OCCT BRepBlend_AppSurface(Funct, First, Last, Tol3d, Tol2d,
    /// TolAngular, Continuity, Degmax, Segmax) (BRepBlend_AppSurface.cxx
    /// L26-95).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        funct: &'a dyn ApproxSweepFunction,
        first: f64,
        last: f64,
        tol3d: f64,
        tol2d: f64,
        tol_angular: f64,
        continuity: crate::geomalgo::gtests_stubs::GeomAbsShape,
        degmax: i32,
        segmax: i32,
    ) -> Self {
        let mut approx = ApproxSweepApproximation::new(funct);
        let nb2d = funct.nb2dcurves();
        let mut nb_pol_sect = 0i32;
        let mut nb_knot_sect = 0i32;
        let mut udeg = 0i32;
        let mut continuity = continuity;

        // (1) Verification de la possibilite de derivation
        if continuity != crate::geomalgo::gtests_stubs::GeomAbsShape::C0 {
            let nb2d = if nb2d == 0 { 1 } else { nb2d };
            funct.section_shape(&mut nb_pol_sect, &mut nb_knot_sect, &mut udeg);
            let n = nb_pol_sect.max(0) as usize;
            let n2d = nb2d.max(0) as usize;
            let mut p = vec![DVec3::ZERO; n];
            let mut p2d = vec![DVec2::ZERO; n2d];
            let mut w = vec![0.0; n];
            // OCCT passes the same buffer twice (V, V) / (W, W); Rust
            // borrows them apart — the duplicate buffers are scratch-only.
            let mut v_a = vec![DVec3::ZERO; n];
            let mut v_b = vec![DVec3::ZERO; n];
            let mut v2d_a = vec![DVec2::ZERO; n2d];
            let mut v2d_b = vec![DVec2::ZERO; n2d];
            let mut w_b = vec![0.0; n];
            let mut w_c = vec![0.0; n];
            let mut ok;
            if continuity == crate::geomalgo::gtests_stubs::GeomAbsShape::C2 {
                ok = funct.d2(
                    first, first, last, &mut p, &mut v_a, &mut v_b, &mut p2d, &mut v2d_a,
                    &mut v2d_b, &mut w, &mut w_b, &mut w_c,
                );
                if !ok {
                    continuity = crate::geomalgo::gtests_stubs::GeomAbsShape::C1;
                }
            }
            if continuity == crate::geomalgo::gtests_stubs::GeomAbsShape::C1 {
                ok = funct.d1(
                    first, first, last, &mut p, &mut v_a, &mut p2d, &mut v2d_a, &mut w, &mut w_b,
                );
                if !ok {
                    continuity = crate::geomalgo::gtests_stubs::GeomAbsShape::C0;
                }
            }
        }

        // (2) Approximation
        // OCCT: approx.Perform(First, Last, Tol3d, Tol3d, Tol2d, TolAngular,
        // continuity, Degmax, Segmax);
        approx.perform(first, last, tol3d, tol3d, tol2d, tol_angular, continuity, degmax, segmax);

        BRepBlendAppSurface { approx }
    }

    /// OCCT BRepBlend_AppSurface(Funct, First, Last, Tol3d, Tol2d,
    /// TolAngular, Continuity) with the OCCT default arguments
    /// Degmax = 11, Segmax = 50 (BRepBlend_AppSurface.hxx L48-57).
    #[allow(clippy::too_many_arguments)]
    pub fn new_defaulted(
        funct: &'a dyn ApproxSweepFunction,
        first: f64,
        last: f64,
        tol3d: f64,
        tol2d: f64,
        tol_angular: f64,
        continuity: crate::geomalgo::gtests_stubs::GeomAbsShape,
    ) -> Self {
        Self::new(funct, first, last, tol3d, tol2d, tol_angular, continuity, 11, 50)
    }
}

impl AppBlendApprox for BRepBlendAppSurface<'_> {
    /// OCCT IsDone() (BRepBlend_AppSurface.lxx L28-31).
    fn is_done(&self) -> bool {
        self.approx.is_done()
    }

    /// OCCT SurfShape(...) (BRepBlend_AppSurface.cxx L97-109).
    fn surf_shape(
        &self,
        u_degree: &mut i32,
        v_degree: &mut i32,
        nb_u_poles: &mut i32,
        nb_v_poles: &mut i32,
        nb_u_knots: &mut i32,
        nb_v_knots: &mut i32,
    ) {
        // OCCT: approx.SurfShape(UDegree, VDegree, NbUPoles, NbVPoles,
        // NbUKnots, NbVKnots); — pending Approx_SweepApproximation; the
        // not-done path reports empty shape.
        *u_degree = 0;
        *v_degree = 0;
        *nb_u_poles = 0;
        *nb_v_poles = 0;
        *nb_u_knots = 0;
        *nb_v_knots = 0;
    }

    /// OCCT Surface(...) (BRepBlend_AppSurface.cxx L111-127).
    #[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
    fn surface(
        &self,
        _t_poles: &mut [Vec<DVec3>],
        _t_weights: &mut [Vec<f64>],
        _t_u_knots: &mut [f64],
        _t_v_knots: &mut [f64],
        _t_u_mults: &mut [i32],
        _t_v_mults: &mut [i32],
    ) {
        // OCCT: approx.Surface(...); — pending Approx_SweepApproximation.
    }

    /// OCCT UDegree() (BRepBlend_AppSurface.lxx L33-36).
    fn u_degree(&self) -> i32 {
        0
    }

    /// OCCT VDegree() (BRepBlend_AppSurface.lxx L38-41).
    fn v_degree(&self) -> i32 {
        0
    }

    /// OCCT SurfPoles() (BRepBlend_AppSurface.lxx L43-46).
    fn surf_poles(&self) -> &[Vec<DVec3>] {
        &[]
    }

    /// OCCT SurfWeights() (BRepBlend_AppSurface.lxx L48-51).
    fn surf_weights(&self) -> &[Vec<f64>] {
        &[]
    }

    /// OCCT SurfUKnots() (BRepBlend_AppSurface.lxx L53-56).
    fn surf_u_knots(&self) -> &[f64] {
        &[]
    }

    /// OCCT SurfVKnots() (BRepBlend_AppSurface.lxx L58-61).
    fn surf_v_knots(&self) -> &[f64] {
        &[]
    }

    /// OCCT SurfUMults() (BRepBlend_AppSurface.lxx L63-66).
    fn surf_u_mults(&self) -> &[i32] {
        &[]
    }

    /// OCCT SurfVMults() (BRepBlend_AppSurface.lxx L68-71).
    fn surf_v_mults(&self) -> &[i32] {
        &[]
    }

    /// OCCT NbCurves2d() (BRepBlend_AppSurface.lxx L73-76).
    fn nb_curves2d(&self) -> i32 {
        0
    }

    /// OCCT Curves2dShape(...) (BRepBlend_AppSurface.cxx L129-133).
    fn curves2d_shape(&self, degree: &mut i32, nb_poles: &mut i32, nb_knots: &mut i32) {
        *degree = 0;
        *nb_poles = 0;
        *nb_knots = 0;
    }

    /// OCCT Curve2d(Index, TPoles, TKnots, TMults)
    /// (BRepBlend_AppSurface.cxx L135-142).
    fn curve2d(&self, _index: i32, _t_poles: &mut [DVec2], _t_knots: &mut [f64], _t_mults: &mut [i32]) {
    }

    /// OCCT Curves2dDegree() (BRepBlend_AppSurface.lxx L78-81).
    fn curves2d_degree(&self) -> i32 {
        0
    }

    /// OCCT Curve2dPoles(Index) (BRepBlend_AppSurface.lxx L83-86).
    fn curve2d_poles(&self, _index: i32) -> &[DVec2] {
        &[]
    }

    /// OCCT Curves2dKnots() (BRepBlend_AppSurface.lxx L88-91).
    fn curves2d_knots(&self) -> &[f64] {
        &[]
    }

    /// OCCT Curves2dMults() (BRepBlend_AppSurface.lxx L93-96).
    fn curves2d_mults(&self) -> &[i32] {
        &[]
    }

    /// OCCT TolReached(Tol3d, Tol2d) (BRepBlend_AppSurface.cxx L144-154).
    fn tol_reached(&self, tol3d: &mut f64, tol2d: &mut f64) {
        // OCCT: Tol3d = approx.MaxErrorOnSurf(); Tol2d = max over the 2d
        // errors; pending Approx_SweepApproximation reports zero errors.
        *tol3d = 0.0;
        *tol2d = 0.0;
    }

    /// OCCT TolCurveOnSurf(Index) (BRepBlend_AppSurface.cxx L140-142).
    fn tol_curve_on_surf(&self, _index: i32) -> f64 {
        0.0
    }
}
