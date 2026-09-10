//! OCCT BRepOffsetAPI_MakePipeShell — 1:1 translation.
//!
//! Source: $OCCT_SRC/src/ModelingAlgorithms/TKOffset/BRepOffsetAPI/
//!         BRepOffsetAPI_MakePipeShell.cxx (L26-267) +
//!         BRepOffsetAPI_MakePipeShell.hxx (L91-395).
//!
//! OCCT inheritance chain (hxx L91): BRepOffsetAPI_MakePipeShell ->
//! BRepPrimAPI_MakeSweep -> BRepBuilderAPI_MakeShape.  Rust has no
//! inheritance: the base-class members (myShape, myGenerated, the Done
//! flag) are kept as plain fields of the struct (the Stage 2e facade
//! precedent).
//!
//! Architecture differences:
//! 1. NCollection_List<TopoDS_Shape> -> Vec<Shape>.
//! 2. The engine member Handle(BRepFill_PipeShell) myPipe (TKBool/BRepFill)
//!    has no rcad translation yet — the BRepFillPipeShell carrier below
//!    keeps the OCCT constructor/method surface with GAP panics (port plan
//!    section 0.6); every facade body around the engine calls is translated
//!    1:1, including the GeomFill_PipeError -> BRepBuilderAPI_PipeError
//!    switch (cxx L157-181).
//! 3. OCCT enums BRepBuilderAPI_PipeError (BRepBuilderAPI_PipeError.hxx),
//!    GeomFill_PipeError (GeomFill_PipeError.hxx),
//!    BRepFill_TypeOfContact (TKBool/BRepFill_TypeOfContact.hxx) and the
//!    BRepBuilderAPI_TransitionMode / BRepFill_TransitionStyle pair are the
//!    local enums below (the transition pair is shared with
//!    brep_offset_api_make_draft.rs).
//! 4. occ::handle<Law_Function> L -> &dyn LawFunction (the landed
//!    geomalgo::law::law_function trait).

use rcad_kernel::math::gp::Ax2;
use rcad_kernel::topo_shape::Shape;

use glam::DVec3;

use crate::geomalgo::law::law_function::LawFunction;

use super::brep_offset_api_make_draft::{
    BRepBuilderAPITransitionMode, BRepFillTransitionStyle,
};

// ---------------------------------------------------------------------------
// OCCT enums (architecture difference #3).
// ---------------------------------------------------------------------------

/// OCCT BRepBuilderAPI_PipeError (BRepBuilderAPI_PipeError.hxx L25-30).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepBuilderAPIPipeError {
    PipeDone,
    PipeNotDone,
    PlaneNotIntersectGuide,
    ImpossibleContact,
}

/// OCCT GeomFill_PipeError (GeomFill_PipeError.hxx L25-34).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeomFillPipeError {
    PipeOk,
    PipeNotDone,
    PlaneNotIntersectGuide,
    ImpossibleContact,
}

/// OCCT BRepFill_TypeOfContact (TKBool/BRepFill_TypeOfContact.hxx L23-27).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BRepFillTypeOfContact {
    NoContact,
    Contact,
    ContactOnBorder,
}

// ---------------------------------------------------------------------------
// GAP carrier (architecture difference #2).
// ---------------------------------------------------------------------------

/// OCCT BRepFill_PipeShell (TKBool/BRepFill, BRepFill_PipeShell.hxx) — the
/// shell-along-a-spine engine of MakePipeShell (architecture difference #2;
/// GAP: no rcad translation yet — the GAP panics are the section 0.6
/// annotation; the constructor and field storage keep the OCCT form).
pub struct BRepFillPipeShell {
    #[allow(dead_code)]
    my_is_done: bool, // OCCT: myIsDone
}

impl BRepFillPipeShell {
    /// OCCT BRepFill_PipeShell::BRepFill_PipeShell(Spine)
    /// (BRepFill_PipeShell.cxx L43) — GAP.
    pub fn new(_the_spine: &Shape) -> Self {
        BRepFillPipeShell { my_is_done: false }
    }

    /// OCCT BRepFill_PipeShell::Set(IsFrenet) — GAP.
    pub fn set(&mut self, _is_frenet: bool) {
        panic!("GAP: BRepFill_PipeShell::Set (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::SetDiscrete() — GAP.
    pub fn set_discrete(&mut self) {
        panic!("GAP: BRepFill_PipeShell::SetDiscrete (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::Set(Axe) — GAP.
    pub fn set_with_axe(&mut self, _axe: &Ax2) {
        panic!("GAP: BRepFill_PipeShell::Set (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::Set(BiNormal) — GAP.
    pub fn set_with_binormal(&mut self, _binormal: &DVec3) {
        panic!("GAP: BRepFill_PipeShell::Set (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::Set(SpineSupport) — GAP.
    pub fn set_with_support(&mut self, _spine_support: &Shape) -> bool {
        panic!("GAP: BRepFill_PipeShell::Set (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::Set(AuxiliarySpine, CurvilinearEquivalence,
    /// KeepContact) — GAP.
    pub fn set_with_auxiliary_spine(
        &mut self,
        _auxiliary_spine: &Shape,
        _curvilinear_equivalence: bool,
        _keep_contact: BRepFillTypeOfContact,
    ) {
        panic!("GAP: BRepFill_PipeShell::Set (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::Add(Profile, WithContact, WithCorrection) —
    /// GAP.
    pub fn add(&mut self, _profile: &Shape, _with_contact: bool, _with_correction: bool) {
        panic!("GAP: BRepFill_PipeShell::Add (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::Add(Profile, Location, WithContact,
    /// WithCorrection) — GAP.
    pub fn add_with_location(
        &mut self,
        _profile: &Shape,
        _location: &Shape,
        _with_contact: bool,
        _with_correction: bool,
    ) {
        panic!("GAP: BRepFill_PipeShell::Add (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::SetLaw(Profile, L, WithContact,
    /// WithCorrection) — GAP.
    pub fn set_law(
        &mut self,
        _profile: &Shape,
        _l: &dyn LawFunction,
        _with_contact: bool,
        _with_correction: bool,
    ) {
        panic!("GAP: BRepFill_PipeShell::SetLaw (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::SetLaw(Profile, L, Location, WithContact,
    /// WithCorrection) — GAP.
    pub fn set_law_with_location(
        &mut self,
        _profile: &Shape,
        _l: &dyn LawFunction,
        _location: &Shape,
        _with_contact: bool,
        _with_correction: bool,
    ) {
        panic!("GAP: BRepFill_PipeShell::SetLaw (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::DeleteProfile(Profile) — GAP.
    pub fn delete_profile(&mut self, _profile: &Shape) {
        panic!("GAP: BRepFill_PipeShell::DeleteProfile (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::IsReady() — GAP.
    pub fn is_ready(&self) -> bool {
        panic!("GAP: BRepFill_PipeShell::IsReady (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::GetStatus() — GAP.
    pub fn get_status(&self) -> GeomFillPipeError {
        panic!("GAP: BRepFill_PipeShell::GetStatus (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::SetTolerance(Tol3d, BoundTol, TolAngular) —
    /// GAP.
    pub fn set_tolerance(&mut self, _tol3d: f64, _bound_tol: f64, _tol_angular: f64) {
        panic!("GAP: BRepFill_PipeShell::SetTolerance (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::SetMaxDegree(NewMaxDegree) — GAP.
    pub fn set_max_degree(&mut self, _new_max_degree: i32) {
        panic!("GAP: BRepFill_PipeShell::SetMaxDegree (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::SetMaxSegments(NewMaxSegments) — GAP.
    pub fn set_max_segments(&mut self, _new_max_segments: i32) {
        panic!("GAP: BRepFill_PipeShell::SetMaxSegments (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::SetForceApproxC1(ForceApproxC1) — GAP.
    pub fn set_force_approx_c1(&mut self, _force_approx_c1: bool) {
        panic!("GAP: BRepFill_PipeShell::SetForceApproxC1 (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::SetTransition(Transition) — GAP.
    pub fn set_transition(&mut self, _transition: BRepFillTransitionStyle) {
        panic!("GAP: BRepFill_PipeShell::SetTransition (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::Simulate(N, R) — GAP.
    pub fn simulate(&mut self, _n: i32, _r: &mut Vec<Shape>) {
        panic!("GAP: BRepFill_PipeShell::Simulate (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::Build() — GAP.
    pub fn build(&mut self) -> bool {
        panic!("GAP: BRepFill_PipeShell::Build (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::MakeSolid() — GAP.
    pub fn make_solid(&mut self) -> bool {
        panic!("GAP: BRepFill_PipeShell::MakeSolid (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::Shape().
    pub fn shape(&self) -> Shape {
        Shape::null()
    }

    /// OCCT BRepFill_PipeShell::FirstShape() — GAP.
    pub fn first_shape(&mut self) -> Shape {
        panic!("GAP: BRepFill_PipeShell::FirstShape (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::LastShape() — GAP.
    pub fn last_shape(&mut self) -> Shape {
        panic!("GAP: BRepFill_PipeShell::LastShape (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::Generated(S, L) — GAP.
    pub fn generated(&mut self, _s: &Shape, _l: &mut Vec<Shape>) {
        panic!("GAP: BRepFill_PipeShell::Generated (TKBool/BRepFill not translated)");
    }

    /// OCCT BRepFill_PipeShell::ErrorOnSurface() — GAP.
    pub fn error_on_surface(&self) -> f64 {
        panic!("GAP: BRepFill_PipeShell::ErrorOnSurface (TKBool/BRepFill not translated)");
    }
}

/// OCCT BRepOffsetAPI_MakePipeShell (hxx L91-395).
pub struct BRepOffsetAPIMakePipeShell {
    // OCCT BRepBuilderAPI base members.
    my_done: bool,            // OCCT BRepBuilderAPI_Command: myDone
    my_shape: Shape,          // OCCT BRepBuilderAPI_MakeShape: myShape
    #[allow(dead_code)]
    my_generated: Vec<Shape>, // OCCT: myGenerated (NCollection_List)
    // OCCT private member (hxx L393): Handle(BRepFill_PipeShell) myPipe.
    my_pipe: BRepFillPipeShell, // OCCT: myPipe
}

impl BRepOffsetAPIMakePipeShell {
    /// OCCT BRepOffsetAPI_MakePipeShell::BRepOffsetAPI_MakePipeShell(Spine)
    /// (cxx L28-34).
    pub fn new(the_spine: &Shape) -> Self {
        // OCCT L30: myPipe = new (BRepFill_PipeShell)(Spine).
        let mut r = BRepOffsetAPIMakePipeShell {
            my_done: false,
            my_shape: Shape::null(),
            my_generated: Vec::new(),
            my_pipe: BRepFillPipeShell::new(the_spine),
        };
        // OCCT L31: SetTolerance().
        r.set_tolerance(1.0e-4, 1.0e-4, 1.0e-2);
        // OCCT L32: SetTransitionMode().
        r.set_transition_mode(BRepBuilderAPITransitionMode::TransitionMode_Transformed);
        // OCCT L33: NotDone().
        r.my_done = false;
        r
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMode(IsFrenet) (cxx L36-39).
    pub fn set_mode(&mut self, is_frenet: bool) {
        self.my_pipe.set(is_frenet);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetDiscreteMode() (cxx L41-44).
    pub fn set_discrete_mode(&mut self) {
        self.my_pipe.set_discrete();
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMode(Axe) (cxx L46-49).
    pub fn set_mode_with_axe(&mut self, axe: &Ax2) {
        self.my_pipe.set_with_axe(axe);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMode(BiNormal) (cxx L51-54).
    pub fn set_mode_with_binormal(&mut self, binormal: &DVec3) {
        self.my_pipe.set_with_binormal(binormal);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMode(SpineSupport)
    /// (cxx L56-59).
    pub fn set_mode_with_support(&mut self, spine_support: &Shape) -> bool {
        self.my_pipe.set_with_support(spine_support)
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMode(AuxiliarySpine,
    /// CurvilinearEquivalence, KeepContact) (cxx L61-67).
    pub fn set_mode_with_auxiliary_spine(
        &mut self,
        auxiliary_spine: &Shape,
        curvilinear_equivalence: bool,
        keep_contact: BRepFillTypeOfContact,
    ) {
        self.my_pipe
            .set_with_auxiliary_spine(auxiliary_spine, curvilinear_equivalence, keep_contact);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Add(Profile, WithContact,
    /// WithCorrection) (cxx L69-74).
    pub fn add(&mut self, profile: &Shape, with_contact: bool, with_correction: bool) {
        self.my_pipe.add(profile, with_contact, with_correction);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Add(Profile, Location, WithContact,
    /// WithCorrection) (cxx L76-82).
    pub fn add_with_location(
        &mut self,
        profile: &Shape,
        location: &Shape,
        with_contact: bool,
        with_correction: bool,
    ) {
        self.my_pipe
            .add_with_location(profile, location, with_contact, with_correction);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetLaw(Profile, L, WithContact,
    /// WithCorrection) (cxx L84-90).
    pub fn set_law(
        &mut self,
        profile: &Shape,
        l: &dyn LawFunction,
        with_contact: bool,
        with_correction: bool,
    ) {
        self.my_pipe.set_law(profile, l, with_contact, with_correction);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetLaw(Profile, L, Location,
    /// WithContact, WithCorrection) (cxx L92-99).
    pub fn set_law_with_location(
        &mut self,
        profile: &Shape,
        l: &dyn LawFunction,
        location: &Shape,
        with_contact: bool,
        with_correction: bool,
    ) {
        self.my_pipe
            .set_law_with_location(profile, l, location, with_contact, with_correction);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Delete(Profile) (cxx L101-105).
    pub fn delete(&mut self, profile: &Shape) {
        self.my_pipe.delete_profile(profile);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::IsReady() (cxx L107-111).
    pub fn is_ready(&self) -> bool {
        self.my_pipe.is_ready()
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::GetStatus() (cxx L113-134).
    pub fn get_status(&self) -> BRepBuilderAPIPipeError {
        // OCCT L116-117: GeomFill_PipeError stat = myPipe->GetStatus().
        let stat = self.my_pipe.get_status();
        // OCCT L118-133: the switch.
        match stat {
            GeomFillPipeError::PipeOk => BRepBuilderAPIPipeError::PipeDone,
            GeomFillPipeError::PlaneNotIntersectGuide => {
                BRepBuilderAPIPipeError::PlaneNotIntersectGuide
            }
            GeomFillPipeError::ImpossibleContact => BRepBuilderAPIPipeError::ImpossibleContact,
            _ => BRepBuilderAPIPipeError::PipeNotDone,
        }
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetTolerance(Tol3d, BoundTol,
    /// TolAngular) (cxx L136-141).
    pub fn set_tolerance(&mut self, tol3d: f64, bound_tol: f64, tol_angular: f64) {
        self.my_pipe.set_tolerance(tol3d, bound_tol, tol_angular);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMaxDegree(NewMaxDegree)
    /// (cxx L143-147).
    pub fn set_max_degree(&mut self, new_max_degree: i32) {
        self.my_pipe.set_max_degree(new_max_degree);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetMaxSegments(NewMaxSegments)
    /// (cxx L149-153).
    pub fn set_max_segments(&mut self, new_max_segments: i32) {
        self.my_pipe.set_max_segments(new_max_segments);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetForceApproxC1(ForceApproxC1)
    /// (cxx L155-165).
    pub fn set_force_approx_c1(&mut self, force_approx_c1: bool) {
        self.my_pipe.set_force_approx_c1(force_approx_c1);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::SetTransitionMode(Mode)
    /// (cxx L167-171).
    pub fn set_transition_mode(&mut self, mode: BRepBuilderAPITransitionMode) {
        // OCCT L169: myPipe->SetTransition((BRepFill_TransitionStyle)Mode) —
        // the enum cast maps TransitionMode_Transformed to
        // TransitionStyle_Transformed (the identical ordinal form).
        let style = match mode {
            BRepBuilderAPITransitionMode::TransitionMode_Transformed => {
                BRepFillTransitionStyle::TransitionStyle_Modified
            }
            BRepBuilderAPITransitionMode::RightCorner => {
                BRepFillTransitionStyle::TransitionStyle_Right
            }
            BRepBuilderAPITransitionMode::RoundCorner => {
                BRepFillTransitionStyle::TransitionStyle_Round
            }
        };
        self.my_pipe.set_transition(style);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Simulate(N, R) (cxx L173-177).
    pub fn simulate(&mut self, n: i32, r: &mut Vec<Shape>) {
        self.my_pipe.simulate(n, r);
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Build(...) (cxx L179-193).
    pub fn build(&mut self) {
        // OCCT L181-182: bool Ok; Ok = myPipe->Build().
        let ok = self.my_pipe.build();
        if ok {
            // OCCT L184-186.
            self.my_shape = self.my_pipe.shape();
            self.my_done = true;
        } else {
            // OCCT L188-191.
            self.my_done = false;
        }
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::MakeSolid() (cxx L195-208).
    pub fn make_solid(&mut self) -> bool {
        // OCCT L197-199: StdFail_NotDone_Raise_if(!IsDone(), ...).
        assert!(
            self.my_done,
            "BRepOffsetAPI_MakePipeShell::MakeSolid"
        );
        // OCCT L200-201: bool Ok; Ok = myPipe->MakeSolid().
        let ok = self.my_pipe.make_solid();
        if ok {
            // OCCT L203: myShape = myPipe->Shape().
            self.my_shape = self.my_pipe.shape();
        }
        // OCCT L207: return Ok.
        ok
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::FirstShape() (cxx L210-213).
    pub fn first_shape(&mut self) -> Shape {
        self.my_pipe.first_shape()
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::LastShape() (cxx L215-218).
    pub fn last_shape(&mut self) -> Shape {
        self.my_pipe.last_shape()
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::Generated(S) (cxx L220-224).
    pub fn generated(&mut self, s: &Shape) -> Vec<Shape> {
        // OCCT L222: myPipe->Generated(S, myGenerated).
        let mut my_generated = std::mem::take(&mut self.my_generated);
        self.my_pipe.generated(s, &mut my_generated);
        self.my_generated = my_generated;
        self.my_generated.clone()
    }

    /// OCCT BRepOffsetAPI_MakePipeShell::ErrorOnSurface() (cxx L226-230).
    pub fn error_on_surface(&self) -> f64 {
        self.my_pipe.error_on_surface()
    }

    /// OCCT BRepBuilderAPI_Command::IsDone() — a PUBLIC member of the OCCT
    /// API (BRepBuilderAPI_Command.hxx: Standard_Boolean IsDone() const);
    /// the rcad field keeps the crate visibility, the accessor restores the
    /// OCCT-public surface (form restoration, no behavior change).
    pub fn is_done(&self) -> bool {
        self.my_done
    }

    /// OCCT BRepBuilderAPI_MakeShape::Shape() — a PUBLIC member of the OCCT
    /// API (BRepBuilderAPI_MakeShape.hxx: Standard_EXPORT const TopoDS_Shape&
    /// Shape() const; raises StdFail_NotDone when not done).  The engine-level
    /// BRepFillPipeShell::shape() accessor above is the OCCT
    /// BRepFill_PipeShell::Shape() member.
    pub fn shape(&self) -> Shape {
        assert!(
            self.my_done,
            "StdFail_NotDone: BRepOffsetAPI_MakePipeShell::Shape()"
        );
        self.my_shape.clone()
    }

    /// Test-world extraction bridge (the FilletResult.brep pattern,
    /// algo_ext::topods_ext::extract_result_brep): flattens the result root
    /// shape into the self-contained BRep the test world consumes
    /// (StepWriter / total_surface_area).  OCCT has no equivalent (the
    /// TopoDS_Shape carries its arena implicitly) — the rcad BRep-pool
    /// architecture difference #4 glue.
    pub fn result_brep(&mut self) -> Option<rcad_kernel::topo::topods::BRep> {
        if !self.my_done || self.my_shape.is_null() {
            return None;
        }
        Some(crate::algo_ext::topods_ext::extract_result_brep(
            &self.my_shape,
            Vec::new(),
        ))
    }
}

