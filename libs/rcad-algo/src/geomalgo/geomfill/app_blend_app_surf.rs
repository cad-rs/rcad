//! OCCT AppBlend_AppSurf (TKGeomAlgo/AppBlend) — 1:1 Rust translation of
//! `AppBlend_AppSurf.gxx` (L1-1134) + `AppBlend_AppSurf.lxx` (L1-170).
//!
//! Template mapping: the OCCT `.gxx` is compiled twice through
//! `#define TheSectionGenerator`:
//! - `GeomFill_AppSurf` (GeomFill_AppSurf.hxx L180-186) over
//!   `GeomFill_SectionGenerator` (rcad [`SectionGenerator`]);
//! - `GeomFill_AppSweep` (GeomFill_AppSweep.hxx L180-186) over
//!   `GeomFill_SweepSectionGenerator` (rcad [`SweepSectionGenerator`]).
//! The template parameter maps to the [`TheSectionGenerator`] trait (the
//! exact OCCT method names); the engine struct carries no generator state —
//! the generator is passed per call — so the two OCCT instantiations
//! collapse into one engine whose Perform entry points are generic over the
//! trait (the ApproxInt_Approx.gxx precedent in `crate::geomalgo::approx_int`).
//! `TheLine` = GeomFill_Line (rcad [`Line`], a plain value — the OCCT
//! `handle` null check of Perform is unrepresentable and kept as a comment).
//!
//! The `UseSmoothing` branch of InternalPerform (gxx L474-541) drives the
//! real `AppDef_Variational` engine
//! ([`crate::geomalgo::app_def_variational::AppDefVariational`], wired in
//! file `app_blend_app_surf_b.rs`); the OCCT
//! `catch (Standard_Failure)` around `Approximate()` maps to the
//! `catch_unwind` convention (hider.rs precedent).

use glam::{DVec2, DVec3};
use rcad_kernel::math::bspl_lib::reparametrize;
use rcad_kernel::math::GeomAbsShape;

use crate::geomalgo::approx_int::ApproxParamType;

use super::line::Line;
use super::section_generator::SectionGenerator;
use super::sweep_section_generator::SweepSectionGenerator;

#[path = "app_blend_app_surf_b.rs"]
mod app_blend_app_surf_b;
#[path = "app_blend_app_surf_c.rs"]
mod app_blend_app_surf_c;

/// OCCT `static bool scal = 1;` (AppBlend_AppSurf.gxx L36) — the file
/// static read by Perform(Lin, F, NbMaxP) (gxx L758); it is never written,
/// so the `if (!scal)` branch is dead in OCCT as well.
const SCAL: bool = true;

/// OCCT Standard_Real RealLast() (Standard_Real.hxx L179-182) == DBL_MAX.
/// OCCT Standard_Real RealFirst() (Standard_Real.hxx L167-170) == -DBL_MAX.
/// — the canonical kernel constants (local duplicates retired).
use rcad_kernel::core::precision::{REAL_FIRST, REAL_LAST};

/// OCCT gp::Resolution() == RealSmall() (gp.hxx L59-60) — the smallest
/// positive normalized double.
const GP_RESOLUTION: f64 = f64::MIN_POSITIVE;

/// OCCT AppBlend_ContextApproxWithNoTgt (AppBlend_Debug.cxx L31-33) — the
/// file static of the debug context switch.
///
/// OCCT AppBlend_AppSurf.gxx L38-39 declares
/// `AppBlend_GetContextSplineApprox` too, but the .gxx body never calls it
/// (only GetContextApproxWithNoTgt at L238/L298/L336/L844); both context
/// pairs are translated from AppBlend_Debug.cxx for completeness.
static APP_BLEND_CONTEXT_SPLINE_APPROX: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);
/// OCCT AppBlend_ContextApproxWithNoTgt (AppBlend_Debug.cxx L37-39).
static APP_BLEND_CONTEXT_APPROX_WITH_NO_TGT: std::sync::atomic::AtomicBool =
    std::sync::atomic::AtomicBool::new(false);

/// OCCT AppBlend_SetContextSplineApprox (AppBlend_Debug.cxx L25-28).
pub fn app_blend_set_context_spline_approx(b: bool) {
    APP_BLEND_CONTEXT_SPLINE_APPROX.store(b, std::sync::atomic::Ordering::Relaxed);
}

/// OCCT AppBlend_GetContextSplineApprox (AppBlend_Debug.cxx L30-33).
#[allow(dead_code)]
pub fn app_blend_get_context_spline_approx() -> bool {
    APP_BLEND_CONTEXT_SPLINE_APPROX.load(std::sync::atomic::Ordering::Relaxed)
}

/// OCCT AppBlend_SetContextApproxWithNoTgt (AppBlend_Debug.cxx L37-40).
pub fn app_blend_set_context_approx_with_no_tgt(b: bool) {
    APP_BLEND_CONTEXT_APPROX_WITH_NO_TGT.store(b, std::sync::atomic::Ordering::Relaxed);
}

/// OCCT AppBlend_GetContextApproxWithNoTgt (AppBlend_Debug.cxx L42-45).
pub fn app_blend_get_context_approx_with_no_tgt() -> bool {
    APP_BLEND_CONTEXT_APPROX_WITH_NO_TGT.load(std::sync::atomic::Ordering::Relaxed)
}

/// OCCT `TheSectionGenerator` — the AppBlend_AppSurf.gxx template parameter
/// (GeomFill_AppSurf.hxx L180 / GeomFill_AppSweep.hxx L180). The trait
/// carries the exact OCCT method names consumed by the .gxx body:
/// GetShape / Knots / Mults / Section / Section / Parameter.
pub trait TheSectionGenerator {
    /// OCCT GetShape(NbPoles, NbKnots, Degree, NbPoles2d).
    fn get_shape(
        &self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    );
    /// OCCT Knots(TKnots).
    fn knots(&self, t_knots: &mut [f64]);
    /// OCCT Mults(TMults).
    fn mults(&self, t_mults: &mut [i32]);
    /// OCCT Section(P, Poles, Poles2d, Weigths).
    fn section(&self, p: i32, poles: &mut [DVec3], poles2d: &mut [DVec2], weigths: &mut [f64]);
    /// OCCT Section(P, Poles, DPoles, Poles2d, DPoles2d, Weigths, DWeigths)
    /// — returns true when the derivatives are computed.
    #[allow(clippy::too_many_arguments)]
    fn section_d1(
        &self,
        p: i32,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool;
    /// OCCT Parameter(P).
    fn parameter(&self, p: i32) -> f64;
}

/// OCCT TheSectionGenerator = GeomFill_SectionGenerator
/// (GeomFill_AppSurf.hxx L180) — pure delegation over the landed
/// geomfill/section_generator.rs translation.
impl TheSectionGenerator for SectionGenerator {
    /// OCCT GeomFill_SectionGenerator::GetShape (cxx L50-60).
    fn get_shape(
        &self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        SectionGenerator::get_shape(self, nb_poles, nb_knots, degree, nb_poles2d);
    }

    /// OCCT GeomFill_SectionGenerator::Knots (cxx L64-67).
    fn knots(&self, t_knots: &mut [f64]) {
        SectionGenerator::knots(self, t_knots);
    }

    /// OCCT GeomFill_SectionGenerator::Mults (cxx L71-74).
    fn mults(&self, t_mults: &mut [i32]) {
        SectionGenerator::mults(self, t_mults);
    }

    /// OCCT GeomFill_SectionGenerator::Section(P, Poles, Poles2d, Weigths)
    /// (cxx L93-102).
    fn section(&self, p: i32, poles: &mut [DVec3], poles2d: &mut [DVec2], weigths: &mut [f64]) {
        SectionGenerator::section(self, p, poles, poles2d, weigths);
    }

    /// OCCT GeomFill_SectionGenerator::Section(P, Poles, DPoles, Poles2d,
    /// DPoles2d, Weigths, DWeigths) (cxx L78-89).
    fn section_d1(
        &self,
        p: i32,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool {
        SectionGenerator::section_d1(self, p, poles, dpoles, poles2d, dpoles2d, weigths, dweigths)
    }

    /// OCCT GeomFill_SectionGenerator::Parameter (cxx L106-109).
    fn parameter(&self, p: i32) -> f64 {
        SectionGenerator::parameter(self, p)
    }
}

/// OCCT TheSectionGenerator = GeomFill_SweepSectionGenerator
/// (GeomFill_AppSweep.hxx L180) — pure delegation over the landed
/// geomfill/sweep_section_generator.rs translation.
impl TheSectionGenerator for SweepSectionGenerator {
    /// OCCT GeomFill_SweepSectionGenerator::GetShape (cxx L376-404).
    fn get_shape(
        &self,
        nb_poles: &mut i32,
        nb_knots: &mut i32,
        degree: &mut i32,
        nb_poles2d: &mut i32,
    ) {
        SweepSectionGenerator::get_shape(self, nb_poles, nb_knots, degree, nb_poles2d);
    }

    /// OCCT GeomFill_SweepSectionGenerator::Knots (cxx L406-428).
    fn knots(&self, t_knots: &mut [f64]) {
        SweepSectionGenerator::knots(self, t_knots);
    }

    /// OCCT GeomFill_SweepSectionGenerator::Mults (cxx L430-450).
    fn mults(&self, t_mults: &mut [i32]) {
        SweepSectionGenerator::mults(self, t_mults);
    }

    /// OCCT GeomFill_SweepSectionGenerator::Section(P, Poles, Poles2d,
    /// Weigths) (cxx L542-667).
    fn section(&self, p: i32, poles: &mut [DVec3], poles2d: &mut [DVec2], weigths: &mut [f64]) {
        SweepSectionGenerator::section(self, p, poles, poles2d, weigths);
    }

    /// OCCT GeomFill_SweepSectionGenerator::Section(P, Poles, DPoles,
    /// Poles2d, DPoles2d, Weigths, DWeigths) (cxx L452-540).
    fn section_d1(
        &self,
        p: i32,
        poles: &mut [DVec3],
        dpoles: &mut [DVec3],
        poles2d: &mut [DVec2],
        dpoles2d: &mut [DVec2],
        weigths: &mut [f64],
        dweigths: &mut [f64],
    ) -> bool {
        SweepSectionGenerator::section_d1(
            self, p, poles, dpoles, poles2d, dpoles2d, weigths, dweigths,
        )
    }

    /// OCCT GeomFill_SweepSectionGenerator::Parameter (cxx L681-698).
    fn parameter(&self, p: i32) -> f64 {
        SweepSectionGenerator::parameter(self, p)
    }
}

/// OCCT AppBlend_AppSurf members (GeomFill_AppSurf.hxx L157-177; the
/// shared .gxx object layout). The NCollection handles map to Option<Vec>
/// (the null-handle state is preserved).
pub struct AppBlendAppSurf {
    /// OCCT bool done.
    pub(super) done: bool,
    /// OCCT int dmin.
    pub(super) dmin: i32,
    /// OCCT int dmax.
    pub(super) dmax: i32,
    /// OCCT double tol3d.
    pub(super) tol3d: f64,
    /// OCCT double tol2d.
    pub(super) tol2d: f64,
    /// OCCT int nbit.
    pub(super) nbit: i32,
    /// OCCT int udeg.
    pub(super) udeg: i32,
    /// OCCT int vdeg.
    pub(super) vdeg: i32,
    /// OCCT bool knownp.
    pub(super) knownp: bool,
    /// OCCT handle(NCollection_HArray2<gp_Pnt>) tabPoles
    /// (1..NbUPoles, 1..NbVPoles) — the [row = U][col = V] rcad form.
    pub(super) tab_poles: Option<Vec<Vec<DVec3>>>,
    /// OCCT handle(NCollection_HArray2<double>) tabWeights.
    pub(super) tab_weights: Option<Vec<Vec<f64>>>,
    /// OCCT handle(NCollection_HArray1<double>) tabUKnots (1..NbUKnots).
    pub(super) tab_u_knots: Option<Vec<f64>>,
    /// OCCT handle(NCollection_HArray1<double>) tabVKnots.
    pub(super) tab_v_knots: Option<Vec<f64>>,
    /// OCCT handle(NCollection_HArray1<int>) tabUMults.
    pub(super) tab_u_mults: Option<Vec<i32>>,
    /// OCCT handle(NCollection_HArray1<int>) tabVMults.
    pub(super) tab_v_mults: Option<Vec<i32>>,
    /// OCCT NCollection_Sequence<handle(NCollection_HArray1<gp_Pnt2d>)>
    /// seqPoles2d.
    pub(super) seq_poles2d: Vec<Vec<DVec2>>,
    /// OCCT double tol3dreached.
    pub(super) tol3dreached: f64,
    /// OCCT double tol2dreached.
    pub(super) tol2dreached: f64,
    /// OCCT Approx_ParametrizationType paramtype.
    pub(super) paramtype: ApproxParamType,
    /// OCCT GeomAbs_Shape continuity.
    pub(super) continuity: GeomAbsShape,
    /// OCCT double critweights[3].
    pub(super) critweights: [f64; 3],
}

impl AppBlendAppSurf {
    /// OCCT AppBlend_AppSurf::AppBlend_AppSurf() (gxx L46-64).
    pub fn new() -> Self {
        AppBlendAppSurf {
            done: false,
            dmin: 0,
            dmax: 0,
            tol3d: 0.0,
            tol2d: 0.0,
            nbit: 0,
            udeg: 0,
            vdeg: 0,
            knownp: false,
            tab_poles: None,
            tab_weights: None,
            tab_u_knots: None,
            tab_v_knots: None,
            tab_u_mults: None,
            tab_v_mults: None,
            seq_poles2d: Vec::new(),
            tol3dreached: 0.0,
            tol2dreached: 0.0,
            paramtype: ApproxParamType::ChordLength,
            continuity: GeomAbsShape::C2,
            critweights: [0.0; 3],
        }
    }

    /// OCCT AppBlend_AppSurf::AppBlend_AppSurf(Degmin, Degmax, Tol3d, Tol2d,
    /// NbIt, KnownParameters) (gxx L68-91).
    pub fn new_with_parameters(
        degmin: i32,
        degmax: i32,
        tol3d: f64,
        tol2d: f64,
        nb_it: i32,
        known_parameters: bool,
    ) -> Self {
        AppBlendAppSurf {
            done: false,
            dmin: degmin,
            dmax: degmax,
            tol3d,
            tol2d,
            nbit: nb_it,
            udeg: 0,
            vdeg: 0,
            knownp: known_parameters,
            tab_poles: None,
            tab_weights: None,
            tab_u_knots: None,
            tab_v_knots: None,
            tab_u_mults: None,
            tab_v_mults: None,
            seq_poles2d: Vec::new(),
            tol3dreached: 0.0,
            tol2dreached: 0.0,
            paramtype: ApproxParamType::ChordLength,
            continuity: GeomAbsShape::C2,
            critweights: [0.0; 3],
        }
    }

    /// OCCT AppBlend_AppSurf::Init (gxx L95-114).
    pub fn init(
        &mut self,
        degmin: i32,
        degmax: i32,
        tol3d: f64,
        tol2d: f64,
        nb_it: i32,
        known_parameters: bool,
    ) {
        self.done = false;
        self.dmin = degmin;
        self.dmax = degmax;
        self.tol3d = tol3d;
        self.tol2d = tol2d;
        self.nbit = nb_it;
        self.knownp = known_parameters;
        self.continuity = GeomAbsShape::C2;
        self.paramtype = ApproxParamType::ChordLength;
        self.critweights[0] = 0.4;
        self.critweights[1] = 0.2;
        self.critweights[2] = 0.4;
    }

    /// OCCT AppBlend_AppSurf::CriteriumWeight (gxx L122-127).
    pub fn criterium_weight(&self) -> (f64, f64, f64) {
        (
            self.critweights[0],
            self.critweights[1],
            self.critweights[2],
        )
    }

    /// OCCT AppBlend_AppSurf::SetCriteriumWeight (gxx L131-138).
    pub fn set_criterium_weight(&mut self, w1: f64, w2: f64, w3: f64) {
        if w1 < 0.0 || w2 < 0.0 || w3 < 0.0 {
            panic!("Standard_DomainError: AppBlend_AppSurf::SetCriteriumWeight");
        }
        self.critweights[0] = w1;
        self.critweights[1] = w2;
        self.critweights[2] = w3;
    }

    /// OCCT AppBlend_AppSurf::SetContinuity (gxx L142-145).
    pub fn set_continuity(&mut self, the_cont: GeomAbsShape) {
        self.continuity = the_cont;
    }

    /// OCCT AppBlend_AppSurf::Continuity (gxx L149-152).
    pub fn continuity(&self) -> GeomAbsShape {
        self.continuity
    }

    /// OCCT AppBlend_AppSurf::SetParType (gxx L156-159).
    pub fn set_par_type(&mut self, par_type: ApproxParamType) {
        self.paramtype = par_type;
    }

    /// OCCT AppBlend_AppSurf::ParType (gxx L163-166).
    pub fn par_type(&self) -> ApproxParamType {
        self.paramtype
    }

    /// OCCT AppBlend_AppSurf::Perform(Lin, F, SpApprox) (gxx L170-176).
    pub fn perform<G: TheSectionGenerator>(
        &mut self,
        lin: &Line,
        f: &mut G,
        sp_approx: bool,
    ) {
        self.internal_perform(lin, f, sp_approx, false);
    }

    /// OCCT AppBlend_AppSurf::PerformSmoothing(Lin, F) (gxx L180-184).
    pub fn perform_smoothing<G: TheSectionGenerator>(&mut self, lin: &Line, f: &mut G) {
        self.internal_perform(lin, f, true, true);
    }

    // -----------------------------------------------------------------
    // The lxx inline queries (AppBlend_AppSurf.lxx L24-170).  The lxx
    // `throw StdFail_NotDone()` / `throw Standard_DomainError()` guards map
    // to the rcad panic convention.
    // -----------------------------------------------------------------

    /// OCCT AppBlend_AppSurf::IsDone (lxx L24-27).
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// OCCT AppBlend_AppSurf::UDegree (lxx L29-36).
    pub fn u_degree(&self) -> i32 {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::UDegree");
        }
        self.udeg
    }

    /// OCCT AppBlend_AppSurf::VDegree (lxx L38-45).
    pub fn v_degree(&self) -> i32 {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::VDegree");
        }
        self.vdeg
    }

    /// OCCT AppBlend_AppSurf::SurfPoles (lxx L47-54).
    pub fn surf_poles(&self) -> &Vec<Vec<DVec3>> {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::SurfPoles");
        }
        self.tab_poles.as_ref().expect("null tabPoles")
    }

    /// OCCT AppBlend_AppSurf::SurfWeights (lxx L56-63).
    pub fn surf_weights(&self) -> &Vec<Vec<f64>> {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::SurfWeights");
        }
        self.tab_weights.as_ref().expect("null tabWeights")
    }

    /// OCCT AppBlend_AppSurf::SurfUKnots (lxx L65-72).
    pub fn surf_u_knots(&self) -> &Vec<f64> {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::SurfUKnots");
        }
        self.tab_u_knots.as_ref().expect("null tabUKnots")
    }

    /// OCCT AppBlend_AppSurf::SurfVKnots (lxx L74-81).
    pub fn surf_v_knots(&self) -> &Vec<f64> {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::SurfVKnots");
        }
        self.tab_v_knots.as_ref().expect("null tabVKnots")
    }

    /// OCCT AppBlend_AppSurf::SurfUMults (lxx L83-90).
    pub fn surf_u_mults(&self) -> &Vec<i32> {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::SurfUMults");
        }
        self.tab_u_mults.as_ref().expect("null tabUMults")
    }

    /// OCCT AppBlend_AppSurf::SurfVMults (lxx L92-99).
    pub fn surf_v_mults(&self) -> &Vec<i32> {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::SurfVMults");
        }
        self.tab_v_mults.as_ref().expect("null tabVMults")
    }

    /// OCCT AppBlend_AppSurf::NbCurves2d (lxx L101-108).
    pub fn nb_curves2d(&self) -> usize {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::NbCurves2d");
        }
        self.seq_poles2d.len()
    }

    /// OCCT AppBlend_AppSurf::Curves2dDegree (lxx L110-121).
    pub fn curves2d_degree(&self) -> i32 {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::Curves2dDegree");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: AppBlend_AppSurf::Curves2dDegree");
        }
        self.vdeg
    }

    /// OCCT AppBlend_AppSurf::Curve2dPoles (lxx L123-134).
    pub fn curve2d_poles(&self, index: i32) -> &Vec<DVec2> {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::Curve2dPoles");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: AppBlend_AppSurf::Curve2dPoles");
        }
        &self.seq_poles2d[(index - 1) as usize]
    }

    /// OCCT AppBlend_AppSurf::Curves2dKnots (lxx L136-147).
    pub fn curves2d_knots(&self) -> &Vec<f64> {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::Curves2dKnots");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: AppBlend_AppSurf::Curves2dKnots");
        }
        self.tab_v_knots.as_ref().expect("null tabVKnots")
    }

    /// OCCT AppBlend_AppSurf::Curves2dMults (lxx L149-160).
    pub fn curves2d_mults(&self) -> &Vec<i32> {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::Curves2dMults");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: AppBlend_AppSurf::Curves2dMults");
        }
        self.tab_v_mults.as_ref().expect("null tabVMults")
    }

    /// OCCT AppBlend_AppSurf::TolReached (lxx L162-170).
    pub fn tol_reached(&self) -> (f64, f64) {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::TolReached");
        }
        (self.tol3dreached, self.tol2dreached)
    }

    // -----------------------------------------------------------------
    // The gxx out-param queries (AppBlend_AppSurf.gxx L1053-1134).
    // -----------------------------------------------------------------

    /// OCCT AppBlend_AppSurf::SurfShape (gxx L1053-1070).  NCollection
    /// Array2 note: ColLength() is the row count and RowLength() the column
    /// count (NCollection_Array2.hxx L198-201).
    #[allow(clippy::too_many_arguments)]
    pub fn surf_shape(
        &self,
        u_degree: &mut i32,
        v_degree: &mut i32,
        nb_u_poles: &mut i32,
        nb_v_poles: &mut i32,
        nb_u_knots: &mut i32,
        nb_v_knots: &mut i32,
    ) {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::SurfShape");
        }
        *u_degree = self.udeg;
        *v_degree = self.vdeg;
        let tab_poles = self.tab_poles.as_ref().expect("null tabPoles");
        *nb_u_poles = tab_poles.len() as i32;
        *nb_v_poles = tab_poles[0].len() as i32;
        *nb_u_knots = self.tab_u_knots.as_ref().expect("null tabUKnots").len() as i32;
        *nb_v_knots = self.tab_v_knots.as_ref().expect("null tabVKnots").len() as i32;
    }

    /// OCCT AppBlend_AppSurf::Surface (gxx L1072-1090) — the array
    /// copy-out form (TPoles = tabPoles->Array2() etc.).
    #[allow(clippy::too_many_arguments)]
    pub fn surface(
        &self,
        t_poles: &mut Vec<Vec<DVec3>>,
        t_weights: &mut Vec<Vec<f64>>,
        t_u_knots: &mut Vec<f64>,
        t_v_knots: &mut Vec<f64>,
        t_u_mults: &mut Vec<i32>,
        t_v_mults: &mut Vec<i32>,
    ) {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::Surface");
        }
        *t_poles = self.tab_poles.as_ref().expect("null tabPoles").clone();
        *t_weights = self.tab_weights.as_ref().expect("null tabWeights").clone();
        *t_u_knots = self.tab_u_knots.as_ref().expect("null tabUKnots").clone();
        *t_u_mults = self.tab_u_mults.as_ref().expect("null tabUMults").clone();
        *t_v_knots = self.tab_v_knots.as_ref().expect("null tabVKnots").clone();
        *t_v_mults = self.tab_v_mults.as_ref().expect("null tabVMults").clone();
    }

    /// OCCT AppBlend_AppSurf::Curves2dShape (gxx L1094-1107).
    pub fn curves2d_shape(&self, degree: &mut i32, nb_poles: &mut i32, nb_knots: &mut i32) {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::Curves2dShape");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: AppBlend_AppSurf::Curves2dShape");
        }
        *degree = self.vdeg;
        *nb_poles = self.tab_poles.as_ref().expect("null tabPoles").len() as i32;
        *nb_knots = self.tab_v_knots.as_ref().expect("null tabVKnots").len() as i32;
    }

    /// OCCT AppBlend_AppSurf::Curve2d (gxx L1111-1127).
    pub fn curve_2d(
        &self,
        index: i32,
        t_poles: &mut Vec<DVec2>,
        t_knots: &mut Vec<f64>,
        t_mults: &mut Vec<i32>,
    ) {
        if !self.done {
            panic!("StdFail_NotDone: AppBlend_AppSurf::Curve2d");
        }
        if self.seq_poles2d.is_empty() {
            panic!("Standard_DomainError: AppBlend_AppSurf::Curve2d");
        }
        *t_poles = self.seq_poles2d[(index - 1) as usize].clone();
        *t_knots = self.tab_v_knots.as_ref().expect("null tabVKnots").clone();
        *t_mults = self.tab_v_mults.as_ref().expect("null tabVMults").clone();
    }

    /// OCCT AppBlend_AppSurf::TolCurveOnSurf (gxx L1131-1134) —
    /// "On ne s'embete pas !!".
    pub fn tol_curve_on_surf(&self, _index: i32) -> f64 {
        self.tol3dreached
    }

    /// OCCT BSplCLib::Reparametrize bridge — used by the Perform bodies
    /// (file-local re-export for the child modules).
    pub(super) fn bspl_reparametrize(u1: f64, u2: f64, knots: &mut [f64]) {
        reparametrize(u1, u2, knots);
    }
}

impl Default for AppBlendAppSurf {
    /// OCCT DEFINE_STANDARD_ALLOC default.
    fn default() -> Self {
        Self::new()
    }
}
